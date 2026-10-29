// src/renderer.ts

declare const electronAPI: {
  openFile: () => Promise<{ fileName: string, data: string } | null>;
};
declare const Plotly: any;

interface SpectroPoint {
  wl: number;
  abs: number;
  smoothAbs?: number;
  tangent?: number;
}

// Constantes issues de votre script R
const CSV_SKIP_ROWS = 16;
const SMOOTH_GAUSSIAN_SIGMA = 3;
const INF_POINT_LEFT_LLOW = 353;
const INF_POINT_LEFT_LHIGH = 473;
const INF_POINT_RIGHT_LLOW = 479;
const INF_POINT_RIGHT_LHIGH = 550;
const CHECK_LEFT = 350;
const CHECK_RIGHT = 750;
const HB_LAMBDA = 414; // Dans votre script c'était 415 pour HB, 414 est classique, j'utilise 415 pour correspondre
const HB_LAMBDA_R = 415;
const BL_LAMBDA = 476;
const RESIDUALS_THRESHOLD_WARNING = 1e-4;
const RESIDUALS_THRESHOLD_ERROR = 1e-3;

document.getElementById('btn-load')?.addEventListener('click', async () => {
  const result = await electronAPI.openFile();
  if (result) {
    document.getElementById('file-name')!.innerText = result.fileName;
    processData(result.data, result.fileName);
  }
});

function ksmooth(x: number[], y: number[], bandwidth: number): number[] {
  // Implémentation basique d'un lissage par noyau gaussien (similaire à ksmooth en R)
  const smoothed: number[] = [];
  const sigma = bandwidth / 2.355; // Approximation FWHM vers sigma
  const variance = sigma * sigma;

  for (let i = 0; i < x.length; i++) {
    const targetX = x[i];
    let weightSum = 0;
    let ySum = 0;

    for (let j = 0; j < x.length; j++) {
      const distance = x[j] - targetX;
      // On ignore les points trop éloignés pour l'optimisation (ex: > 3 sigma)
      if (Math.abs(distance) > 3 * sigma) continue;
      
      const weight = Math.exp(-(distance * distance) / (2 * variance));
      weightSum += weight;
      ySum += y[j] * weight;
    }
    smoothed.push(weightSum === 0 ? y[i] : ySum / weightSum);
  }
  return smoothed;
}

function parseTstData(rawData: string): SpectroPoint[] {
  // Séparer par tout type d'espace (espaces multiples, sauts de ligne, etc.)
  const parts = rawData.split(/\s+/);
  
  // Convertir toutes les chaînes en nombres (les textes deviendront NaN, ce qui est attendu)
  const fileValues = parts.map(p => parseFloat(p));

  let readIndex = 0;
  let lastValidTrace: SpectroPoint[] = [];
  let nTraces = 0;

  while (readIndex < fileValues.length - 1) {
    const curValue = fileValues[readIndex];
    const nextValue = fileValues[readIndex + 1];

    // On cherche deux nombres valides consécutifs
    if (!isNaN(curValue) && !isNaN(nextValue)) {
      // Détection du début de la trace (3500 = 350.0 nm)
      if (curValue === 3500 && nextValue !== 7500) {
        let traceValid = true;
        const currentTrace: SpectroPoint[] = [];

        // Boucle sur les longueurs d'onde attendues : 350 à 750 par pas de 2
        for (let traceIndex = 350; traceIndex <= 750; traceIndex += 2) {
          if (readIndex + 1 >= fileValues.length) {
            console.warn("Le fichier se termine au milieu d'une trace.");
            traceValid = false;
            break;
          }

          const wl = fileValues[readIndex] / 10; // ex: 3500 -> 350
          const t = fileValues[readIndex + 1];   // transmittance

          if (isNaN(wl) || isNaN(t)) {
            console.warn(`Valeur NaN trouvée à la position: ${readIndex}, attendu: ${traceIndex}`);
            traceValid = false;
            break;
          }

          if (wl !== traceIndex) {
            console.warn(`Lambda incohérent à la position: ${readIndex}, attendu: ${traceIndex}, reçu: ${wl}`);
            traceValid = false;
            break;
          }

          // Calcul de l'absorbance à partir de la transmittance
          const abs = Math.log10(1 / t);
          currentTrace.push({ wl, abs });
          
          readIndex += 2;
        }

        // Contrôle qualité de la trace
        if (traceValid && currentTrace.length === 201) {
          lastValidTrace = currentTrace;
          nTraces++;
        }
        // Le `for` a déjà avancé le readIndex, on passe à l'itération suivante du `while`
        continue;
      }
    }
    readIndex++;
  }

  if (nTraces === 0 || lastValidTrace.length === 0) {
    throw new Error("Données invalides : impossible de charger la trace depuis le fichier .tst");
  }

  // Comme dans le script R, on renvoie la dernière trace valide trouvée dans le fichier
  return lastValidTrace;
}

function processData(rawData: string, fileName: string) {
  let spectrum: SpectroPoint[] = [];

  try {
    if (fileName.toLowerCase().endsWith('.tst')) {
      // ---------------------------------------------------------
      // Utilisation du parseur .tst
      // ---------------------------------------------------------
      spectrum = parseTstData(rawData);
      
    } else if (fileName.toLowerCase().endsWith('.txt') || fileName.toLowerCase().endsWith('.csv')) {
      // Parseur texte générique (CSV/TXT 2 colonnes)
      const lines = rawData.split('\n');
      for (const line of lines) {
        const parts = line.split(/[,;\t\s]+/);
        if (parts.length >= 2) {
          const wl = parseFloat(parts[0]);
          const abs = parseFloat(parts[1]);
          if (!isNaN(wl) && !isNaN(abs) && wl >= 350 && wl <= 750) {
            spectrum.push({ wl, abs });
          }
        }
      }
    } else {
      alert("Format de fichier non pris en charge par ce parseur (utilisez .tst, .csv ou .txt).");
      return;
    }
  } catch (error: any) {
    alert(error.message);
    return;
  }

  // Trier par longueur d'onde pour sécuriser la suite
  spectrum.sort((a, b) => a.wl - b.wl);

  // 1. Lissage des données
  const wls = spectrum.map(p => p.wl);
  const rawAbs = spectrum.map(p => p.abs);
  const smoothAbsArray = ksmooth(wls, rawAbs, SMOOTH_GAUSSIAN_SIGMA);
  
  for (let i = 0; i < spectrum.length; i++) {
    spectrum[i].smoothAbs = smoothAbsArray[i];
  }

  // 2. Recherche de la meilleure tangente (get_auto_best_tangent_parameters)
  let minResiduals = 999;
  let bestLeftWl = NaN;
  let bestRightWl = NaN;
  let bestTangentValues: number[] = [];

  const leftCandidates = spectrum.filter(p => p.wl >= INF_POINT_LEFT_LLOW && p.wl <= INF_POINT_LEFT_LHIGH);
  const rightCandidates = spectrum.filter(p => p.wl >= INF_POINT_RIGHT_LLOW && p.wl <= INF_POINT_RIGHT_LHIGH);

  for (const leftPt of leftCandidates) {
    for (const rightPt of rightCandidates) {
      if (rightPt.wl <= leftPt.wl) continue;

      const a = (rightPt.smoothAbs! - leftPt.smoothAbs!) / (rightPt.wl - leftPt.wl);
      const b = rightPt.smoothAbs! - a * rightPt.wl;

      let sumResiduals = 0;
      const currentTangentValues: number[] = [];

      for (const pt of spectrum) {
        const tgValue = a * pt.wl + b;
        currentTangentValues.push(tgValue);
        
        if (pt.wl >= CHECK_LEFT && pt.wl <= CHECK_RIGHT) {
            const residual = pt.smoothAbs! - tgValue;
            if (residual < 0) {
                sumResiduals += (-residual);
            }
        }
      }

      if (sumResiduals < minResiduals) {
        minResiduals = sumResiduals;
        bestLeftWl = leftPt.wl;
        bestRightWl = rightPt.wl;
        bestTangentValues = currentTangentValues;
      }
    }
  }

  // Assignation de la tangente
  for (let i = 0; i < spectrum.length; i++) {
      spectrum[i].tangent = bestTangentValues[i];
  }

  // 3. Calculs quantitatifs (compute_quantitative_results)
  // Recherche des points spécifiques
  const getPtByWl = (targetWl: number) => {
      // Trouve le point le plus proche si l'exact n'existe pas
      return spectrum.reduce((prev, curr) => 
          Math.abs(curr.wl - targetWl) < Math.abs(prev.wl - targetWl) ? curr : prev
      );
  }

  const hbPt = getPtByWl(HB_LAMBDA_R);
  const blPt = getPtByWl(BL_LAMBDA);

  const hbTraceDo = hbPt.smoothAbs!;
  const blTraceDo = blPt.smoothAbs!;
  const hbTangentDo = hbPt.tangent!;
  const blTangentDo = blPt.tangent!;

  const hbDeltaDo = hbTraceDo - hbTangentDo;
  const blDeltaDo = blTraceDo - blTangentDo;

  // 4. Évaluation des conclusions
  const hbConcl = hbDeltaDo >= 0.1 ? ">0.1" : (hbDeltaDo >= 0.02 ? ">0.02" : "<0.02");
  const blConcl = blDeltaDo >= 0.007 ? ">0.007" : "<0.007"; // Note: J'ai corrigé le text R '>0.07' de votre code en '>0.007' pour rester cohérent

  let comConcl = "?";
  if (blConcl === "<0.007") {
    if (hbConcl === ">0.1") {
      comConcl = "C1 : Présence d’hémoglobine en quantité trop importante pouvant masquer la présence de bilirubine. Hémorragie sous-arachnoïdienne non exclue.";
    } else if (hbConcl === ">0.02") {
      comConcl = "C3 : Présence d’hémoglobine en faible quantité mais absence de bilirubine. Pas d’argument en faveur d’une hémorragie sous-arachnoïdienne ou compatible avec une HSA débutante ou une PL traumatique.";
    } else {
      comConcl = "C2 : Examen spectrophotométrique négatif. Absence d’hémoglobine et de bilirubine. Pas d’argument en faveur d’une hémorragie sous-arachnoïdienne.";
    }
  } else {
    if (hbConcl === ">0.1" || hbConcl === ">0.02") {
      comConcl = "C4 : Présence de bilirubine et d’hémoglobine. Résultat compatible avec une hémorragie sous-arachnoïdienne. (Attention: non spécifique, autres étiologies possibles).";
    } else {
      comConcl = "Delta DO bilirubine positif! Regarder la protéinorachie et la bilirubinémie plasmatique pour conclure.";
    }
  }

  // 5. Mise à jour de l'UI
  document.getElementById('val-noa')!.innerText = hbDeltaDo.toFixed(4);
  document.getElementById('val-nba')!.innerText = blDeltaDo.toFixed(4);

  const conclusionEl = document.getElementById('val-conclusion')!;
  conclusionEl.innerText = comConcl;
  
  if (blDeltaDo > 0.007) {
    conclusionEl.style.backgroundColor = '#fee2e2';
    conclusionEl.style.color = '#991b1b';
  } else {
    conclusionEl.style.backgroundColor = '#dcfce7';
    conclusionEl.style.color = '#166534';
  }

  // 6. Tracé Plotly
  const traceRaw = {
    x: spectrum.map(p => p.wl),
    y: spectrum.map(p => p.abs),
    mode: 'lines',
    name: 'Brut',
    line: { color: 'lightgray' }
  };

  const traceSmooth = {
    x: spectrum.map(p => p.wl),
    y: spectrum.map(p => p.smoothAbs),
    mode: 'lines',
    name: 'Lissé',
    line: { color: 'black' }
  };

  const traceTangent = {
    x: spectrum.map(p => p.wl),
    y: spectrum.map(p => p.tangent),
    mode: 'lines',
    name: 'Tangente',
    line: { color: 'red' }
  };

  // Points de la tangente
  const traceTangentPoints = {
      x: [bestLeftWl, bestRightWl],
      y: [getPtByWl(bestLeftWl).smoothAbs, getPtByWl(bestRightWl).smoothAbs],
      mode: 'markers',
      name: 'Points Tangente',
      marker: { color: 'red', size: 8 }
  }

  const layout = {
    title: 'Spectre d\'absorption',
    xaxis: { title: 'Longueur d\'onde (nm)', range: [350, 750] },
    yaxis: { title: 'Absorbance', autorange: true }, // Vous pouvez adapter ylim_max comme dans plot_spectrum
    margin: { t: 40, l: 50, r: 20, b: 40 },
    shapes: [
        // Ligne verticale pour HB
        {
            type: 'line',
            x0: HB_LAMBDA_R, y0: 0, x1: HB_LAMBDA_R, y1: getPtByWl(HB_LAMBDA_R).abs,
            line: { color: 'black', dash: 'dot', width: 2 }
        },
        // Ligne verticale pour BL
        {
            type: 'line',
            x0: BL_LAMBDA, y0: 0, x1: BL_LAMBDA, y1: getPtByWl(BL_LAMBDA).abs,
            line: { color: 'black', dash: 'dot', width: 2 }
        }
    ]
  };

  Plotly.newPlot('plotly-div', [traceRaw, traceSmooth, traceTangent, traceTangentPoints], layout, { responsive: true });
}
