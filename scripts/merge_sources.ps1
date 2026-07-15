# Merge .c files by module; strip duplicate #include "..." from parts after the first.
param(
    [string]$SrcDir = "$PSScriptRoot\..\src"
)

function Merge-CModule {
    param(
        [string]$OutName,
        [string[]]$Parts
    )
    $outPath = Join-Path $SrcDir $OutName
    $sb = [System.Text.StringBuilder]::new()
    [void]$sb.AppendLine("// SPDX-License-Identifier: BSD-3-Clause")
    [void]$sb.AppendLine("/* merged: $($Parts -join ' + ') */")
    [void]$sb.AppendLine()

    $first = $true
    foreach ($part in $Parts) {
        $path = Join-Path $SrcDir $part
        if (-not (Test-Path $path)) { throw "Missing $path" }
        $lines = Get-Content $path -Encoding UTF8
        $i = 0
        while ($i -lt $lines.Count) {
            $l = $lines[$i]
            if ($l -match '^// SPDX') { $i++; continue }
            if ($l.Trim() -eq '/*') {
                $i++
                while ($i -lt $lines.Count -and $lines[$i] -notmatch '^\s*\*/') { $i++ }
                if ($i -lt $lines.Count) { $i++ }
                continue
            }
            break
        }
        if (-not $first) {
            [void]$sb.AppendLine("/* --- $part --- */")
            [void]$sb.AppendLine()
        }
        for (; $i -lt $lines.Count; $i++) {
            $l = $lines[$i]
            if (-not $first -and $l -match '^\s*#\s*include\s+"') { continue }
            [void]$sb.AppendLine($l)
        }
        $first = $false
        [void]$sb.AppendLine()
    }
    [System.IO.File]::WriteAllText($outPath, $sb.ToString(), [System.Text.UTF8Encoding]::new($false))
    Write-Host "Wrote $outPath"
}

$SrcDir = (Resolve-Path $SrcDir).Path

Merge-CModule "fft.c" @("hsfft.c", "real.c")
Merge-CModule "ts.c" @("conv.c", "filter.c", "talg.c")
Merge-CModule "cstats.c" @("erfunc.c", "dist.c", "pdist.c", "stats.c", "polyroot.c")
Merge-CModule "linalg.c" @("matrix.c", "lls.c", "regression.c", "nls.c")
Merge-CModule "optim.c" @("lnsrchmp.c", "newtonmin.c", "conjgrad.c", "secant.c", "brent.c", "neldermead.c", "optimc.c")
Merge-CModule "fit.c" @("emle.c")
