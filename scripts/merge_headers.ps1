# Merge headers into layer *.h (fft.h, cstats.h, ...); one include guard per output file.
param([string]$SrcDir = "$PSScriptRoot\..\src")

function Merge-Headers {
    param([string]$OutName, [string[]]$Parts, [string]$Guard)
    $outPath = Join-Path $SrcDir $OutName
    $sb = [System.Text.StringBuilder]::new()
    [void]$sb.AppendLine("// SPDX-License-Identifier: BSD-3-Clause")
    [void]$sb.AppendLine("/* Module header: merged from $($Parts -join ', ') */")
    [void]$sb.AppendLine("#ifndef $Guard")
    [void]$sb.AppendLine("#define $Guard")
    [void]$sb.AppendLine()
    [void]$sb.AppendLine("#ifdef __cplusplus")
    [void]$sb.AppendLine('extern "C" {')
    [void]$sb.AppendLine("#endif")
    [void]$sb.AppendLine()

    $first = $true
    foreach ($part in $Parts) {
        $path = Join-Path $SrcDir $part
        $text = Get-Content $path -Raw -Encoding UTF8
        # strip license and include guard blocks
        $text = $text -replace '(?s)// SPDX.*?\n', ''
        $text = $text -replace '(?s)/\*.*?\*/\s*', '', 1
        $text = $text -replace '(?s)#ifndef\s+\w+_H_\s*#define\s+\w+_H_\s*', ''
        $text = $text -replace '(?s)#ifdef __cplusplus\s*extern "C" \{\s*#endif\s*', ''
        $text = $text -replace '(?s)#ifdef __cplusplus\s*\}\s*#endif\s*', ''
        $text = $text -replace '#endif\s*/\*\s*\w+_H_\s*\*/\s*', ''
        $text = $text -replace '#endif\s*/\*.*?\*/\s*', ''
        if (-not $first) {
            $text = $text -replace '(?m)^\s*#\s*include\s+"[^"]+"\s*\r?\n', ''
        }
        $text = $text.Trim()
        [void]$sb.AppendLine("/* --- $part --- */")
        [void]$sb.AppendLine($text)
        [void]$sb.AppendLine()
        $first = $false
    }
    [void]$sb.AppendLine("#ifdef __cplusplus")
    [void]$sb.AppendLine("}")
    [void]$sb.AppendLine("#endif")
    [void]$sb.AppendLine()
    [void]$sb.AppendLine("#endif /* $Guard */")
    [void]$sb.AppendLine()
    [System.IO.File]::WriteAllText($outPath, $sb.ToString(), [System.Text.UTF8Encoding]::new($false))
    Write-Host "Wrote $outPath"
}

$SrcDir = (Resolve-Path $SrcDir).Path
Merge-Headers "fft.h" @("hsfft.h", "real.h") "FFT_H_"
Merge-Headers "cstats.h" @("erfunc.h", "dist.h", "pdist.h", "stats.h", "polyroot.h") "CSTATS_H_"
Merge-Headers "linalg.h" @("matrix.h", "lls.h", "regression.h", "nls.h") "LINALG_H_"
Merge-Headers "optim.h" @("lnsrchmp.h", "newtonmin.h", "conjgrad.h", "secant.h", "brent.h", "neldermead.h", "optimc.h") "OPTIM_H_"
Merge-Headers "ts.h" @("conv.h", "filter.h", "talg.h") "TS_H_"
Copy-Item (Join-Path $SrcDir "emle.h") (Join-Path $SrcDir "fit.h") -Force
Write-Host "Wrote fit.h (from emle.h)"
