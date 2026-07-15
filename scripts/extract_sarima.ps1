# Regenerate slim optim.c / linalg.c / cstats.c from ctsa-master full sources.
$ErrorActionPreference = "Stop"
$Src = Join-Path $PSScriptRoot "..\src"
$Master = Join-Path $PSScriptRoot "..\..\ctsa-master\src"

function Get-Functions($text) {
    $lines = $text -split "`n", -1
    $funcs = @{}
    $i = 0
    while ($i -lt $lines.Count) {
        $line = $lines[$i]
        if ($line -match '^(static\s+)?(?:inline\s+)?(?:const\s+)?(?:unsigned\s+|signed\s+)?(?:struct\s+\w+\s*\*?\s*|(?:double|int|void|reg_object|nls_object|opt_object)\s*\*?\s*)*(\w+)\s*\(') {
            $name = $Matches[2]
            if ($name -in @('if','while','for','switch','return')) { $i++; continue }
            $sig = $line
            $j = $i
            while ($sig -notmatch '\)\s*\{' -and $j + 1 -lt $lines.Count) {
                $j++
                $sig += " " + $lines[$j]
            }
            if ($sig -notmatch '\)\s*\{') { $i++; continue }
            $depth = ($sig.ToCharArray() | Where-Object { $_ -eq '{' }).Count - ($sig.ToCharArray() | Where-Object { $_ -eq '}' }).Count
            $k = $j + 1
            while ($k -lt $lines.Count -and $depth -gt 0) {
                $depth += ($lines[$k].ToCharArray() | Where-Object { $_ -eq '{' }).Count - ($lines[$k].ToCharArray() | Where-Object { $_ -eq '}' }).Count
                $k++
            }
            if (-not $funcs.ContainsKey($name)) { $funcs[$name] = @($i, $k) }
            $i = $k
            continue
        }
        $i++
    }
    return @{ Lines = $lines; Funcs = $funcs }
}

function Extract-File($masterName, $keep, $outName, $header) {
    $path = Join-Path $Master $masterName
    if (-not (Test-Path $path)) { throw "Missing upstream $path" }
    $text = [IO.File]::ReadAllText($path)
    $parsed = Get-Functions $text
    $lines = $parsed.Lines
    $funcs = $parsed.Funcs
    foreach ($k in $keep) {
        if (-not $funcs.ContainsKey($k)) { throw "Missing $k in $masterName" }
    }
    $out = New-Object System.Collections.Generic.List[string]
    [void]$out.Add("// SPDX-License-Identifier: BSD-3-Clause")
    [void]$out.Add("/*")
    [void]$out.Add(" * $outName - $header")
    [void]$out.Add(" */")
    [void]$out.Add("")
    if ($outName -eq 'optim.c') { [void]$out.Add('#include "optim.h"') }
    elseif ($outName -eq 'linalg.c') { [void]$out.Add('#include "linalg.h"') }
    else { [void]$out.Add('#include "cstats.h"'); [void]$out.Add('#include "linalg.h"') }
    [void]$out.Add("")
    foreach ($name in $keep) {
        $r = $funcs[$name]
        for ($idx = $r[0]; $idx -lt $r[1]; $idx++) { [void]$out.Add($lines[$idx].TrimEnd("`r")) }
        [void]$out.Add("")
    }
    $outPath = Join-Path $Src $outName
    [IO.File]::WriteAllText($outPath, ($out -join "`n") + "`n", [Text.UTF8Encoding]::new($false))
    Write-Host "Wrote $outName"
}

$optimKeep = @(
    'grad_fd','grad_calc','linsolve_lower','hessian_fd','lnsrch','stopcheck',
    'jrotate','qrupdate','inithess_lower','bfgs_factored','bfgs_min','fminunc'
)
$linalgKeep = @(
    'macheps','signx','imax','imin','l2norm','array_max_abs',
    'scale','madd','msub','nmult','stranspose',
    'pludecomp','ludecomp','linsolve','minverse'
)
$cstatsKeep = @(
    'mean','mcon','CMOD','CDIVID','RESCALE','CAUCHY','ERREV','POLYEV','NEXTH',
    'CALCT','VRSHFT','FXSHFT','NOSHFT','cpoly','polyroot'
)

# Upstream: newtonmin.c+optimc.c -> optim; matrix+lls+regression -> linalg; dist+erfunc+polyroot -> cstats
Extract-File 'newtonmin.c' $optimKeep 'optim.c' 'SARIMA BFGS + fminunc + hessian_fd'
Extract-File 'matrix.c' $linalgKeep 'linalg.c' 'SARIMA matrix ops + LU + inverse'

# mean lives in stats.c upstream; polyroot in polyroot.c
$cstatsPath = Join-Path $Src 'cstats.c'
$meanText = [IO.File]::ReadAllText((Join-Path $Master 'stats.c'))
$polyText = [IO.File]::ReadAllText((Join-Path $Master 'polyroot.c'))
$meanF = Get-Functions $meanText
$polyF = Get-Functions $polyText
foreach ($k in $cstatsKeep) {
    if ($k -eq 'mean') { if (-not $meanF.Funcs.ContainsKey('mean')) { throw 'mean missing' } }
    elseif (-not $polyF.Funcs.ContainsKey($k)) { throw "Missing $k in polyroot.c" }
}
# Build cstats.c manually from stats mean + polyroot helpers
$out = New-Object System.Collections.Generic.List[string]
[void]$out.Add("// SPDX-License-Identifier: BSD-3-Clause")
[void]$out.Add("/*")
[void]$out.Add(" * cstats.c - SARIMA stats (mean, polyroot)")
[void]$out.Add(" */")
[void]$out.Add("")
[void]$out.Add('#include "cstats.h"')
[void]$out.Add('#include "linalg.h"')
[void]$out.Add("")
$r = $meanF.Funcs['mean']
for ($idx = $r[0]; $idx -lt $r[1]; $idx++) { [void]$out.Add($meanF.Lines[$idx].TrimEnd("`r")) }
[void]$out.Add("")
foreach ($name in $cstatsKeep) {
    if ($name -eq 'mean') { continue }
    $r = $polyF.Funcs[$name]
    for ($idx = $r[0]; $idx -lt $r[1]; $idx++) { [void]$out.Add($polyF.Lines[$idx].TrimEnd("`r")) }
    [void]$out.Add("")
}
[IO.File]::WriteAllText($cstatsPath, ($out -join "`n") + "`n", [Text.UTF8Encoding]::new($false))
Write-Host 'Wrote cstats.c'

# SARIMA-specific patches on linalg.c
$lin = [IO.File]::ReadAllText((Join-Path $Src 'linalg.c'))
$lin = $lin -replace '(?s)void mmult\(double\* A, double\* B, double\* C,int m,int n, int p\) \{.*?\n\}', @'
void mmult(double* A, double* B, double* C, int m, int n, int p)
{
	nmult(A, B, C, m, n, p);
}
'@
if ($lin -notmatch 'void mtranspose') {
    $lin += @'

void mtranspose(double *sig, int rows, int cols, double *col)
{
	stranspose(sig, rows, cols, col);
}
'@
}
[IO.File]::WriteAllText((Join-Path $Src 'linalg.c'), $lin, [Text.UTF8Encoding]::new($false))
Write-Host 'Patched linalg.c (mmult, mtranspose)'
