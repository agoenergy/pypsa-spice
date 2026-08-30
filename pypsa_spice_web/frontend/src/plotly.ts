interface PlotlyApi {
  react: (element: HTMLElement, data: unknown[], layout: object, config: object) => void;
  purge: (element: HTMLElement) => void;
  Plots: { resize: (element: HTMLElement) => void };
}

declare global {
  interface Window {
    Plotly?: PlotlyApi;
  }
}

let plotlyPromise: Promise<PlotlyApi> | null = null;

export function loadPlotly(): Promise<PlotlyApi> {
  if (window.Plotly) return Promise.resolve(window.Plotly);
  if (plotlyPromise) return plotlyPromise;

  plotlyPromise = new Promise<PlotlyApi>((resolve, reject) => {
    const script = document.createElement("script");
    script.src = "/vendor/plotly.min.js";
    script.async = true;
    script.dataset.plotlyLoader = "true";
    script.addEventListener("load", () => {
      if (window.Plotly) resolve(window.Plotly);
      else reject(new Error("The Plotly script loaded without exposing its chart API."));
    }, { once: true });
    script.addEventListener("error", () => reject(new Error("The Plotly chart library could not be loaded.")), { once: true });
    document.head.append(script);
  }).catch((reason) => {
    plotlyPromise = null;
    throw reason;
  });

  return plotlyPromise;
}

