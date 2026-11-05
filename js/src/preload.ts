import { contextBridge, ipcRenderer } from 'electron';

contextBridge.exposeInMainWorld('electronAPI', {
  openFile: () => ipcRenderer.invoke('dialog:openFile'),
  exportPdf: (defaultFilename: string) => ipcRenderer.invoke('dialog:exportPdf', defaultFilename)
});
