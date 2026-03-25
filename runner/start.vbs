' Pipeline Runner — silent launcher
' Double-click or use desktop shortcut to start

Dim fso, oShell, repoRoot

Set fso    = CreateObject("Scripting.FileSystemObject")
Set oShell = CreateObject("WScript.Shell")

' Derive repo root from this script's location (runner\ is one level inside repo)
repoRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))

oShell.CurrentDirectory = repoRoot

' Start server hidden (window style 0 = invisible)
oShell.Run "py runner\app.py", 0, False

' Wait for server to be ready
WScript.Sleep 3000

' Open browser
oShell.Run "http://127.0.0.1:8501", 1, False
