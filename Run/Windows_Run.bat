@echo off
TITLE Metadoon
:: Esconde a janela preta se possível, mas mantém logs visíveis por enquanto
echo Starting Metadoon...

:: 1. Volta para a raiz
cd ..\..

:: 2. Ativa o ambiente e roda o Python
call conda run -n metadoon python metadoon.py --no-console

if %errorlevel% neq 0 (
    echo [ERROR] Metadoon crashed or environment not found.
    echo Did you run the Installer first?
    pause
)