#!/bin/bash

# USAGE
#    chmod +x InitProject.sh
#    ./InitProject.sh

# Create Python virtual environment
python3 -m venv venv

# Activate the virtual environment
source ./venv/bin/activate

# Create project folder structure
mkdir -p src data tests scripts docs examples .vscode

# Initialize Git repository
git init

# Create necessary empty files
touch README.md requirements.txt .gitignore .vscode/settings.json

# Write contents to .gitignore
cat << EOF > .gitignore
# Byte-compiled / optimized / DLL files
__pycache__/
*.py[cod]
*$py.class

# C extensions
*.so

# Distribution / packaging
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
pip-wheel-metadata/
share/python-wheels/
*.egg-info/
.installed.cfg
*.egg

# PyInstaller
*.manifest
*.spec

# Installer logs
pip-log.txt
pip-delete-this-directory.txt

# Unit test / coverage reports
htmlcov/
.tox/
.nox/
.coverage
.coverage.*
.cache
nosetests.xml
coverage.xml
*.cover
.hypothesis/
.pytest_cache/

# Translations
*.mo
*.pot

# Django stuff
*.log
local_settings.py
db.sqlite3
db.sqlite3-journal

# Flask stuff
instance/
.webassets-cache

# Scrapy stuff
.scrapy

# Sphinx documentation
docs/_build/

# PyBuilder
target/

# Jupyter Notebook
.ipynb_checkpoints

# pyenv
.python-version

# celery beat schedule file
celerybeat-schedule

# SageMath parsed files
*.sage.py

# Environments
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# Spyder project settings
.spyderproject
.spyproject

# Rope project settings
.ropeproject

# mkdocs documentation
/site

# mypy
.mypy_cache/
.dmypy.json
dmypy.json

# Pyre type checker
.pyre/
EOF

# Write VSCode settings
cat << EOF > .vscode/settings.json
{
    "python.pythonPath": "./venv/bin/python3"
}
EOF

# Completion message
echo "Python project folder has been successfully created."
