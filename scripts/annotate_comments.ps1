# ASCII-only driver: strip comments, add UTF-8 headers from scripts/headers/*.hdr
param([string]$Root = (Join-Path $PSScriptRoot ".."))
$ErrorActionPreference = "Stop"
$Utf8 = New-Object System.Text.UTF8Encoding $false
$ScriptDir = $PSScriptRoot

function Remove-CComments([string]$src) {
    $sb = New-Object System.Text.StringBuilder
    $i = 0; $n = $src.Length; $state = "code"
    while ($i -lt $n) {
        $c = $src[$i]
        $c2 = if ($i + 1 -lt $n) { $src[$i + 1] } else { [char]0 }
        switch ($state) {
            "code" {
                if ($c -eq '"' -and ($i -eq 0 -or $src[$i - 1] -ne '\')) { $state = "str"; [void]$sb.Append($c); $i++; continue }
                if ($c -eq "'" -and ($i -eq 0 -or $src[$i - 1] -ne '\')) { $state = "chr"; [void]$sb.Append($c); $i++; continue }
                if ($c -eq '/' -and $c2 -eq '/') { $state = "line"; $i += 2; continue }
                if ($c -eq '/' -and $c2 -eq '*') { $state = "block"; $i += 2; continue }
                [void]$sb.Append($c); $i++
            }
            "str" { [void]$sb.Append($c); if ($c -eq '"' -and $src[$i - 1] -ne '\') { $state = "code" }; $i++ }
            "chr" { [void]$sb.Append($c); if ($c -eq "'" -and $src[$i - 1] -ne '\') { $state = "code" }; $i++ }
            "line" { if ($c -eq "`n") { $state = "code"; [void]$sb.Append($c) }; $i++ }
            "block" { if ($c -eq '*' -and $c2 -eq '/') { $state = "code"; $i += 2 } else { $i++ } }
        }
    }
    return $sb.ToString()
}

$FuncDescMap = @{}
$descPath = Join-Path $ScriptDir "func_desc.txt"
if (Test-Path $descPath) {
    foreach ($line in [IO.File]::ReadAllLines($descPath, $Utf8)) {
        $line = $line.Trim()
        if ($line -eq "" -or $line.StartsWith("#")) { continue }
        $p = $line.IndexOf("|")
        if ($p -gt 0) { $FuncDescMap[$line.Substring(0, $p)] = $line.Substring($p + 1) }
    }
}

$FuncRules = New-Object System.Collections.Generic.List[object]
$rulesPath = Join-Path $ScriptDir "func_rules.txt"
if (Test-Path $rulesPath) {
    foreach ($line in [IO.File]::ReadAllLines($rulesPath, $Utf8)) {
        $line = $line.Trim()
        if ($line -eq "" -or $line.StartsWith("#")) { continue }
        $p = $line.IndexOf([char]9)
        if ($p -lt 0) { $p = $line.IndexOf("|") }
        if ($p -gt 0) {
            $FuncRules.Add([pscustomobject]@{
                Pattern = $line.Substring(0, $p)
                Prefix  = $line.Substring($p + 1)
            })
        }
    }
}

$defPath = Join-Path $ScriptDir "func_default.txt"
$DefaultPrefix = if (Test-Path $defPath) { [IO.File]::ReadAllText($defPath, $Utf8).Trim() } else { "numeric" }

function Escape-CommentText([string]$s) {
    if ($null -eq $s) { return "" }
    return $s.Replace('*/', '* /')
}

function Get-FuncDescZh([string]$name) {
    if ($FuncDescMap.ContainsKey($name)) { return $FuncDescMap[$name] }
    $n = $name.ToLower()
    foreach ($r in $FuncRules) {
        if ($n -match $r.Pattern) { return $r.Prefix + ": " + $name }
    }
    return $DefaultPrefix + ": " + $name
}

function Test-FuncStartLine([string]$line, [ref]$name) {
    if ($line -match '^\t') { return $false }
    if ($line -match '^    ') { return $false }
    if ($line -match '^\s*#\s*(include|define|if|elif|else|endif|pragma)') { return $false }
    if ($line -match '^\s*(struct|typedef|enum)\b') { return $false }
    if ($line -match '^\s*(if|while|for|switch|else|do|return|case|default)\b') { return $false }
    if ($line -notmatch '\(') { return $false }
    if ($line -match '=\s*\{') { return $false }
    if ($line -match '^\s*(?:static\s+)?(?:const\s+)?(?:unsigned\s+)?(?:signed\s+)?(?:[\w][\w\s\*]*?\b)(\w+)\s*\(') {
        $name.Value = $Matches[1]
        if ($name.Value -in @('if','while','for','switch','else','sizeof','define')) { return $false }
        return $true
    }
    return $false
}

function Find-Functions([string[]]$lines) {
    $funcs = New-Object System.Collections.Generic.List[object]
    $i = 0
    while ($i -lt $lines.Count) {
        $fname = $null
        if (-not (Test-FuncStartLine $lines[$i] ([ref]$fname))) { $i++; continue }
        $sig = $lines[$i].TrimEnd()
        $j = $i
        while ($sig -notmatch '\)\s*\{' -and $j + 1 -lt $lines.Count) {
            $nxt = $lines[$j + 1].Trim()
            if ($nxt -match '^\{') { break }
            if ($nxt -match '^(if|while|for|switch|else)\b') { break }
            $j++; $sig = "$sig $nxt"
        }
        if (($sig -match '\)\s*\{') -or ($j + 1 -lt $lines.Count -and $lines[$j + 1] -match '^\s*\{')) {
            $funcs.Add([pscustomobject]@{ Index = $i; Name = $fname })
            $i = $j + 1; continue
        }
        $i++
    }
    return $funcs
}

function Annotate-CFile([string]$path) {
    $base = [IO.Path]::GetFileName($path)
    $hdrPath = Join-Path $ScriptDir "headers\$base.hdr"
    if (-not (Test-Path $hdrPath)) { Write-Warning "No header for $base"; return }
    $fileHdr = [IO.File]::ReadAllText($hdrPath, $Utf8).TrimEnd()
    $lines = @((Remove-CComments ([IO.File]::ReadAllText($path, $Utf8))) -split "`r?`n")
    $funcs = Find-Functions $lines
    $start = 0
    for ($s = 0; $s -lt $lines.Count; $s++) { if ($lines[$s] -match '^\s*#\s*include') { $start = $s; break } }
    $out = [System.Collections.Generic.List[string]]::new()
    [void]$out.Add("// SPDX-License-Identifier: BSD-3-Clause")
    foreach ($hl in ($fileHdr -split "`r?`n")) { [void]$out.Add($hl) }
    [void]$out.Add("")
    for ($k = $start; $k -lt $lines.Count; $k++) {
        $hit = $funcs | Where-Object { $_.Index -eq $k } | Select-Object -First 1
        if ($hit) {
            $zh = Escape-CommentText (Get-FuncDescZh $hit.Name)
            [void]$out.Add("/*")
            [void]$out.Add(" * " + $hit.Name + " - " + $zh)
            [void]$out.Add(" */")
        }
        [void]$out.Add($lines[$k])
    }
    [IO.File]::WriteAllText($path, (($out -join "`n").TrimEnd() + "`n"), $Utf8)
    Write-Host ("OK " + $base + " n=" + $funcs.Count)
}

@(
    'src\sarima.c','src\ts.c','src\cstats.c',
    'src\linalg.c','src\optim.c','src\fit.c','test\sarimatest.c'
) | ForEach-Object {
    $p = Join-Path $Root $_
    if (Test-Path $p) { Annotate-CFile $p }
}
