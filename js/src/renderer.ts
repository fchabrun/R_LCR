// Déclaration globale pour TypeScript (fournie par preload.js et le CDN Plotly)
declare const electronAPI: {
  openFile: () => Promise<{ fileName: string, data: string } | null>;
};
declare const Plotly: any;

interface SpectroPoint {
  wl: number; // Wavelength
  abs: number; // Absorbance
}

document.getElementById('btn-load')?.addEventListener('click', async () => {
  const result = await electronAPI.openFile();
  if (result) {
    document.getElementById('file-name')!.innerText = result.fileName;
    processData(result.data);
  }
});

function processData(rawData: string) {
  // 1. Parsing simpliste (à adapter selon le format exact de votre spectrophotomètre)
  const lines = rawData.split('\n');
  const spectrum: SpectroPoint[] = [];
  
  for (const line of lines) {
    const [wlStr, absStr] = line.split(/[,;\t]/);
    const wl = parseFloat(wlStr);
    const abs = parseFloat(absStr);
    if (!isNaN(wl) && !isNaN(abs)) {
      spectrum.push({ wl, abs });
    }
  }

  // Trier par longueur d'onde croissante
  spectrum.sort((a, b) => a.wl - b.wl);

  // 2. Calculs (À calquer sur la logique exacte de votre repo R_LCR)
  // Simulation de la ligne de base entre 350 nm et 600 nm
  const pt350 = spectrum.find(p => p.wl >= 350) || spectrum[0];
  const pt600 = spectrum.find(p => p.wl >= 600) || spectrum[spectrum.length - 1];
  
  // Fonction de la droite de base : y = mx + b
  const m = (pt600.abs - pt350.abs) / (pt600.wl - pt350.wl);
  const b = pt350.abs - m * pt350.wl;
  const calcBaseline = (wl: number) => m * wl + b;

  // Récupération des absorbances aux pics clés
  const pt414 = spectrum.find(p => Math.abs(p.wl - 414) <= 1) || { wl: 414, abs: 0 };
  const pt476 = spectrum.find(p => Math.abs(p.wl - 476) <= 1) || { wl: 476, abs: 0 };

  // Calcul NOA et NBA
  const noa = pt414.abs - calcBaseline(414);
  const nba = pt476.abs - calcBaseline(476);

  // 3. Mise à jour de l'UI
  document.getElementById('val-noa')!.innerText = noa.toFixed(4);
  document.getElementById('val-nba')!.innerText = nba.toFixed(4);

  const conclusionEl = document.getElementById('val-conclusion')!;
  if (nba > 0.007) { // Seuil usuel, à adapter
    conclusionEl.style.backgroundColor = '#fee2e2'; // Rouge clair
    conclusionEl.style.color = '#991b1b';
    conclusionEl.innerText = "Profil en faveur d'une HSA (NBA positif)";
  } else {
    conclusionEl.style.backgroundColor = '#dcfce7'; // Vert clair
    conclusionEl.style.color = '#166534';
    conclusionEl.innerText = "Pas de xanthochromie significative";
  }

  // 4. Tracé Plotly
  const traceSpectrum = {
    x: spectrum.map(p => p.wl),
    y: spectrum.map(p => p.abs),
    mode: 'lines',
    name: 'Spectre LCR',
    line: { color: '#2563eb' }
  };

  const traceBaseline = {
    x: [pt350.wl, pt600.wl],
    y: [pt350.abs, pt600.abs],
    mode: 'lines',
    name: 'Ligne de base',
    line: { color: '#9ca3af', dash: 'dash' }
  };

  const layout = {
    title: 'Spectre d\'absorption',
    xaxis: { title: 'Longueur d\'onde (nm)', range: [300, 650] },
    yaxis: { title: 'Absorbance' },
    margin: { t: 40, l: 50, r: 20, b: 40 }
  };

  Plotly.newPlot('plotly-div', [traceSpectrum, traceBaseline], layout, { responsive: true });
}
