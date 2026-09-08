param (
    [Parameter(Mandatory=$false)]
    [string]$Preset = "1"
)

$presets = @{
    "1" = "1_modern_tech.yml"
    "2" = "2_minimalist_gitbook.yml"
    "3" = "3_sci_lab_dark.yml"
    "4" = "4_classic_readthedocs.yml"
    "5" = "5_bootswatch_flatly.yml"
}

if (-not $presets.ContainsKey($Preset)) {
    Write-Host "Available presets:" -ForegroundColor Yellow
    Write-Host "  1 - Modern Tech (Top Tabs, Inter font, Teal accent, Developer style)"
    Write-Host "  2 - Minimalist Clean (GitBook/VitePress style, Integrated TOC, Outfit font)"
    Write-Host "  3 - Sci-Fi / Neuro Lab (Dark Mode default, Cyan/Emerald glow, Lexend font)"
    Write-Host "  4 - Classic ReadTheDocs (Sphinx style)"
    Write-Host "  5 - Bootswatch Flatly (Bootstrap clean style)"
    Write-Host ""
    Write-Host "Usage: .\switch_theme.ps1 <1-5>" -ForegroundColor Cyan
    exit 1
}

$target = Join-Path "theme_presets" $presets[$Preset]
if (Test-Path $target) {
    Copy-Item -Path $target -Destination "mkdocs.yml" -Force
    Write-Host "Applied Preset $Preset ($($presets[$Preset])) -> mkdocs.yml" -ForegroundColor Green
    Write-Host "If 'mkdocs serve' is running, your browser at http://127.0.0.1:8000 will reload automatically." -ForegroundColor Cyan
} else {
    Write-Host "Error: Could not find preset file $target" -ForegroundColor Red
}
