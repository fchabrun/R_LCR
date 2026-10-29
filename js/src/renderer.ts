declare const electronAPI: {
  openFile: () => Promise<{ fileName: string, data: string } | null>;
};
declare const Plotly: any;

interface SpectroPoint { wl: number; abs: number; smoothAbs?: number; tangent?: number; }

const CSV_SKIP_ROWS = 16;
const SMOOTH_GAUSSIAN_SIGMA = 3;
const INF_POINT_LEFT_LLOW = 353;
const INF_POINT_LEFT_LHIGH = 473;
const INF_POINT_RIGHT_LLOW = 479;
const INF_POINT_RIGHT_LHIGH = 550;
const CHECK_LEFT = 350;
const CHECK_RIGHT = 750;
const HB_LAMBDA = 414; // HB_LAMBDA_R
const BL_LAMBDA = 476;
const WARNING_LEVEL_OK = 0;
const WARNING_LEVEL_WARN = 1;
const WARNING_LEVEL_STOP = 2;
const RESIDUALS_THRESHOLD_WARNING = 1e-4;
const RESIDUALS_THRESHOLD_ERROR = 1e-3;

// État global pour permettre le recalcul lors du changement auto/manuel sans recharger le fichier
let currentSpectrum: SpectroPoint[] = [];
let currentRawDataExists = false;

// ---- Événements UI ----
document.getElementById('btn-load')?.addEventListener('click', async () => {
  const result = await electronAPI.openFile();
  if (result) {
    document.getElementById('file-name')!.innerText = result.fileName;
    loadAndProcessFile(result.data, result.fileName);
  }
});

document.getElementsByName('calc-mode').forEach(radio => {
  radio.addEventListener('change', (e) => {
    const isManual = (e.target as HTMLInputElement).value === 'manual';
    document.getElementById('manual-inputs')!.style.display = isManual ? 'flex' : 'none';
    if (currentRawDataExists) updateAnalysis();
  });
});

document.getElementById('manual-left')?.addEventListener('input', () => { if (currentRawDataExists) updateAnalysis(); });
document.getElementById('manual-right')?.addEventListener('input', () => { if (currentRawDataExists) updateAnalysis(); });


// ---- Fonctions Utilitaires Mathématiques ----

function ksmooth(x: number[], y: number[], bandwidth: number): number[] {
  const smoothed: number[] = [];
  const sigma = bandwidth / 2.355;
  const variance = sigma * sigma;

  for (let i = 0; i < x.length; i++) {
    const targetX = x[i];
    let weightSum = 0;
    let ySum = 0;

    for (let j = 0; j < x.length; j++) {
      const distance = x[j] - targetX;
      if (Math.abs(distance) > 3 * sigma) continue;
      const weight = Math.exp(-(distance * distance) / (2 * variance));
      weightSum += weight;
      ySum += y[j] * weight;
    }
    smoothed.push(weightSum === 0 ? y[i] : ySum / weightSum);
  }
  return smoothed;
}

// Équivalent de approx() en R pour avoir des pas de 1 nm constants
function interpolateTo1nm(rawData: SpectroPoint[]): SpectroPoint[] {
  if (rawData.length === 0) return [];
  const minWl = Math.ceil(rawData[0].wl);
  const maxWl = Math.floor(rawData[rawData.length - 1].wl);
  const interpolated: SpectroPoint[] = [];

  for (let wl = minWl; wl <= maxWl; wl++) {
    const exact = rawData.find(p => p.wl === wl);
    if (exact) {
      interpolated.push({ wl, abs: exact.abs });
    } else {
      const left = [...rawData].reverse().find(p => p.wl < wl);
      const right = rawData.find(p => p.wl > wl);
      if (left && right) {
        const ratio = (wl - left.wl) / (right.wl - left.wl);
        const abs = left.abs + ratio * (right.abs - left.abs);
        interpolated.push({ wl, abs });
      }
    }
  }
  return interpolated;
}

// Parseur TST
function parseTstData(rawData: string): SpectroPoint[] {
  const parts = rawData.split(/\s+/);
  const fileValues = parts.map(p => parseFloat(p));
  let readIndex = 0;
  let lastValidTrace: SpectroPoint[] = [];

  while (readIndex < fileValues.length - 1) {
    const curValue = fileValues[readIndex];
    const nextValue = fileValues[readIndex + 1];
    if (!isNaN(curValue) && !isNaN(nextValue)) {
      if (curValue === 3500 && nextValue !== 7500) {
        let traceValid = true;
        const currentTrace: SpectroPoint[] = [];
        for (let traceIndex = 350; traceIndex <= 750; traceIndex += 2) {
          if (readIndex + 1 >= fileValues.length) { traceValid = false; break; }
          const wl = fileValues[readIndex] / 10;
          const t = fileValues[readIndex + 1];
          if (isNaN(wl) || isNaN(t) || wl !== traceIndex) { traceValid = false; break; }
          
          currentTrace.push({ wl, abs: Math.log10(1 / t) });
          readIndex += 2;
        }
        if (traceValid && currentTrace.length === 201) lastValidTrace = currentTrace;
        continue;
      }
    }
    readIndex++;
  }
  if (lastValidTrace.length === 0) throw new Error("Données .tst invalides.");
  return lastValidTrace;
}


// ---- Logique Principale ----

function loadAndProcessFile(rawData: string, fileName: string) {
  let rawSpectrum: SpectroPoint[] = [];

  try {
    if (fileName.toLowerCase().endsWith('.tst')) {
      rawSpectrum = parseTstData(rawData);
    } else if (fileName.toLowerCase().endsWith('.txt') || fileName.toLowerCase().endsWith('.csv')) {
      const lines = rawData.split('\n');
      for (const line of lines) {
        const parts = line.split(/[,;\t\s]+/);
        if (parts.length >= 2) {
          const wl = parseFloat(parts[0]);
          const abs = parseFloat(parts[1]);
          if (!isNaN(wl) && !isNaN(abs) && wl >= 350 && wl <= 750) rawSpectrum.push({ wl, abs });
        }
      }
    } else {
      alert("Format non pris en charge.");
      return;
    }
  } catch (error: any) {
    alert(error.message); return;
  }

  rawSpectrum.sort((a, b) => a.wl - b.wl);
  
  // Interpolation à 1nm (approx dans app.R)
  currentSpectrum = interpolateTo1nm(rawSpectrum);

  // Lissage
  const wls = currentSpectrum.map(p => p.wl);
  const rawAbs = currentSpectrum.map(p => p.abs);
  const smoothAbsArray = ksmooth(wls, rawAbs, SMOOTH_GAUSSIAN_SIGMA);
  currentSpectrum.forEach((pt, i) => pt.smoothAbs = smoothAbsArray[i]);

  currentRawDataExists = true;
  updateAnalysis(); // Lance le calcul et l'affichage
}

function updateAnalysis() {
  const mode = (document.querySelector('input[name="calc-mode"]:checked') as HTMLInputElement).value;
  let leftWl = NaN, rightWl = NaN;

  // -- 1. Calcul de la Tangente --
  if (mode === 'auto') {
    let minResiduals = Infinity;
    const leftCandidates = currentSpectrum.filter(p => p.wl >= INF_POINT_LEFT_LLOW && p.wl <= INF_POINT_LEFT_LHIGH);
    const rightCandidates = currentSpectrum.filter(p => p.wl >= INF_POINT_RIGHT_LLOW && p.wl <= INF_POINT_RIGHT_LHIGH);

    for (const leftPt of leftCandidates) {
      for (const rightPt of rightCandidates) {
        if (rightPt.wl <= leftPt.wl) continue;
        const a = (rightPt.smoothAbs! - leftPt.smoothAbs!) / (rightPt.wl - leftPt.wl);
        const b = rightPt.smoothAbs! - a * rightPt.wl;
        let sumResiduals = 0;
        
        for (const pt of currentSpectrum) {
          const tgValue = a * pt.wl + b;
          if (pt.wl >= CHECK_LEFT && pt.wl <= CHECK_RIGHT) {
            const residual = pt.smoothAbs! - tgValue;
            if (residual < 0) sumResiduals += (-residual);
          }
        }

        if (sumResiduals < minResiduals) {
          minResiduals = sumResiduals;
          leftWl = leftPt.wl;
          rightWl = rightPt.wl;
        }
      }
    }
  } else {
    // Mode manuel
    leftWl = parseFloat((document.getElementById('manual-left') as HTMLInputElement).value);
    rightWl = parseFloat((document.getElementById('manual-right') as HTMLInputElement).value);
  }

  // Application de la droite (y = ax + b)
  const ptLeft = currentSpectrum.find(p => p.wl === leftWl);
  const ptRight = currentSpectrum.find(p => p.wl === rightWl);
  
  if (ptLeft && ptRight) {
    const a = (ptRight.smoothAbs! - ptLeft.smoothAbs!) / (ptRight.wl - ptLeft.wl);
    const b = ptRight.smoothAbs! - a * ptRight.wl;
    currentSpectrum.forEach(pt => pt.tangent = a * pt.wl + b);
  }

  // -- 2. Évaluation Clinique et Qualité (compute_quantitative_results) --
  let warningLevel = WARNING_LEVEL_OK;
  const warningItems: string[] = [];

  if (currentSpectrum.some(p => p.abs < 0)) {
    warningItems.push("Valeurs d'absorbance < 0 !");
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_WARN);
  }
  if (currentSpectrum.some(p => p.abs < 0 && p.wl < 550)) {
    warningItems.push("Valeurs d'absorbance < 0 entre 350 et 550nm !");
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_STOP);
  }
  if (mode !== "auto") {
    warningItems.push("Méthode manuelle choisie: vérifier la cohérence des résultats !");
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_WARN);
  }

  // lambda max où la tangente est positive
  const positiveTangentsWl = currentSpectrum.filter(p => p.tangent! > 0).map(p => p.wl);
  const maxTangentCrossWl = positiveTangentsWl.length > 0 ? Math.max(...positiveTangentsWl) : 0;
  
  if (maxTangentCrossWl < BL_LAMBDA) {
    warningItems.push(`La tangente passe par zéro avant ${BL_LAMBDA}nm !`);
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_WARN);
  }

  // Calcul strict des résidus (tangente qui coupe la courbe)
  const restrictedSpectrum = currentSpectrum.filter(p => p.wl >= INF_POINT_LEFT_LLOW && p.wl <= INF_POINT_RIGHT_LHIGH);
  let sumTotalResiduals = 0;
  for (const pt of restrictedSpectrum) {
    const res = pt.smoothAbs! - pt.tangent!;
    if (res < 0) sumTotalResiduals += (-res); // Si smooth - tg < 0 => tg est au dessus
  }

  if (sumTotalResiduals > RESIDUALS_THRESHOLD_ERROR) {
    warningItems.push(mode === "auto" 
      ? `Impossible de trouver une tangente sans couper la courbe entre ${INF_POINT_LEFT_LLOW} et ${INF_POINT_RIGHT_LHIGH}nm !` 
      : `Les positions choisies sont invalides: la tangente passe au-dessus de la courbe entre ${INF_POINT_LEFT_LLOW} et ${INF_POINT_RIGHT_LHIGH}nm !`);
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_STOP);
  } else if (sumTotalResiduals > RESIDUALS_THRESHOLD_WARNING) {
    warningItems.push(mode === "auto" 
      ? `Impossible de trouver une tangente sans couper la courbe entre ${INF_POINT_LEFT_LLOW} et ${INF_POINT_RIGHT_LHIGH}nm !` 
      : `Les positions choisies sont douteuses: la tangente passe au-dessus de la courbe entre ${INF_POINT_LEFT_LLOW} et ${INF_POINT_RIGHT_LHIGH}nm !`);
    warningLevel = Math.max(warningLevel, WARNING_LEVEL_WARN);
  }

  // Points et Deltas
  const hbPt = currentSpectrum.find(p => p.wl === HB_LAMBDA) || currentSpectrum[0];
  const blPt = currentSpectrum.find(p => p.wl === BL_LAMBDA) || currentSpectrum[0];

  const hbDeltaDo = hbPt.smoothAbs! - hbPt.tangent!;
  const blDeltaDo = blPt.smoothAbs! - blPt.tangent!;

  const hbConcl = hbDeltaDo >= 0.1 ? ">0.1" : (hbDeltaDo >= 0.02 ? ">0.02" : "<0.02");
  const blConcl = blDeltaDo >= 0.007 ? ">0.007" : "<0.007";

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
      comConcl = "C4 : Présence de bilirubine et d’hémoglobine. Résultat compatible avec une hémorragie sous-arachnoïdienne. (Attention: non spécifique).";
    } else {
      comConcl = "Delta DO bilirubine positif! Regarder la protéinorachie et la bilirubinémie plasmatique pour conclure.";
    }
  }

  // -- 3. Mise à jour Interface UI --
  const warningsBox = document.getElementById('warnings-box')!;
  const warningsList = document.getElementById('warnings-list')!;
  
  if (warningItems.length > 0) {
    warningsBox.style.display = 'block';
    warningsBox.className = `warnings-box level-${warningLevel}`;
    warningsList.innerHTML = warningItems.map(w => `<li>${w}</li>`).join('');
  } else {
    warningsBox.style.display = 'none';
  }

  document.getElementById('val-noa')!.innerText = hbDeltaDo.toFixed(4);
  document.getElementById('val-nba')!.innerText = blDeltaDo.toFixed(4);

  const conclusionEl = document.getElementById('val-conclusion')!;
  conclusionEl.innerText = comConcl;
  if (blDeltaDo > 0.007) {
    conclusionEl.style.backgroundColor = '#fee2e2'; conclusionEl.style.color = '#991b1b';
  } else {
    conclusionEl.style.backgroundColor = '#dcfce7'; conclusionEl.style.color = '#166534';
  }

  // -- 4. Tracé Plotly --
  const traceRaw = { x: currentSpectrum.map(p => p.wl), y: currentSpectrum.map(p => p.abs), mode: 'lines', name: 'Brut', line: { color: 'lightgray' } };
  const traceSmooth = { x: currentSpectrum.map(p => p.wl), y: currentSpectrum.map(p => p.smoothAbs), mode: 'lines', name: 'Lissé', line: { color: 'black' } };
  const traceTangent = { x: currentSpectrum.map(p => p.wl), y: currentSpectrum.map(p => p.tangent), mode: 'lines', name: 'Tangente', line: { color: 'red' } };
  
  const traceTangentPoints = {
      x: [leftWl, rightWl],
      y: [ptLeft?.smoothAbs, ptRight?.smoothAbs],
      mode: 'markers', name: 'Points Tangente', marker: { color: 'red', size: 8 }
  };

  const layout = {
    title: 'Spectre d\'absorption',
    xaxis: { title: 'Longueur d\'onde (nm)', range: [350, 750] },
    yaxis: { title: 'Absorbance' },
    margin: { t: 40, l: 50, r: 20, b: 40 },
    shapes: [
        { type: 'line', x0: HB_LAMBDA, y0: 0, x1: HB_LAMBDA, y1: hbPt?.abs || 0.1, line: { color: 'black', dash: 'dot', width: 2 } },
        { type: 'line', x0: BL_LAMBDA, y0: 0, x1: BL_LAMBDA, y1: blPt?.abs || 0.1, line: { color: 'black', dash: 'dot', width: 2 } }
    ]
  };

  Plotly.newPlot('plotly-div', [traceRaw, traceSmooth, traceTangent, traceTangentPoints], layout, { responsive: true });
}
