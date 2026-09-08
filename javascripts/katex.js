function renderKatex() {
  if (typeof renderMathInElement === 'function') {
    renderMathInElement(document.body, {
      delimiters: [
        { left: "$$", right: "$$", display: true },
        { left: "$", right: "$", display: false },
        { left: "\\(", right: "\\)", display: false },
        { left: "\\[", right: "\\]", display: true }
      ],
      throwOnError: false
    });
  }
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", renderKatex);
} else {
  renderKatex();
}

if (typeof document$ !== "undefined") {
  document$.subscribe(renderKatex);
}
