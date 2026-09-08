function initMermaid() {
  if (typeof mermaid === "undefined") {
    console.warn("Mermaid.js script not loaded yet.");
    return;
  }

  mermaid.initialize({
    startOnLoad: false,
    theme: "base",
    themeVariables: {
      fontFamily: '"Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
      fontSize: "14px",
      primaryColor: "#f8fafc",
      primaryTextColor: "#0f172a",
      primaryBorderColor: "#94a3b8",
      lineColor: "#475569",
      secondaryColor: "#f1f5f9",
      tertiaryColor: "#ffffff",
      clusterBkg: "#f8fafc",
      clusterBorder: "#cbd5e1",
      edgeLabelBackground: "#ffffff",
      nodeBorder: "#94a3b8",
      mainBkg: "#f8fafc",
      nodeTextColor: "#0f172a"
    },
    securityLevel: "loose",
    flowchart: {
      useMaxWidth: true,
      htmlLabels: true,
      curve: "basis",
      nodeSpacing: 35,
      rankSpacing: 35,
      padding: 12
    }
  });

  // Find all mermaid code blocks (either <pre class="mermaid"><code> or <div class="mermaid">)
  var elements = document.querySelectorAll("pre.mermaid, div.mermaid, .language-mermaid");
  if (elements.length === 0) return;

  elements.forEach(function (el) {
    // If it's a <pre><code> block, unwrap to a clean div
    if (el.tagName.toLowerCase() === "pre" || el.querySelector("code")) {
      var code = el.querySelector("code") ? el.querySelector("code").textContent : el.textContent;
      var newDiv = document.createElement("div");
      newDiv.className = "mermaid";
      newDiv.textContent = code;
      el.parentNode.replaceChild(newDiv, el);
    }
  });

  mermaid.run({
    querySelector: ".mermaid"
  });
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", initMermaid);
} else {
  initMermaid();
}

if (typeof document$ !== "undefined") {
  document$.subscribe(initMermaid);
}


