#!/bin/bash

# GO Term Analysis - Dependency Installation Script
# This script installs all required Python packages for the GO analysis system

echo "GO Term Analysis - Dependency Installation"
echo "=========================================="

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    echo "Please install Python 3.10+ and try again"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
echo "Found Python version: $PYTHON_VERSION"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "Error: pip3 is not installed or not in PATH"
    echo "Please install pip3 and try again"
    exit 1
fi

echo "Found pip3: $(pip3 --version)"

# Check if virtual environment should be created
read -p "Do you want to create a virtual environment? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Creating virtual environment..."
    
    # Check if venv module is available
    if python3 -c "import venv" &> /dev/null; then
        python3 -m venv go_analysis_env
        echo "Virtual environment created: go_analysis_env"
        
        # Activate virtual environment
        source go_analysis_env/bin/activate
        echo "Virtual environment activated"
        
        # Upgrade pip
        pip install --upgrade pip
    else
        echo "Warning: venv module not available, using system Python"
    fi
else
    echo "Using system Python installation"
fi

echo ""
echo "Installing required packages..."
echo "=============================="

# Install packages from requirements file
if [ -f "requirement/requirements_go_analysis.txt" ]; then
    echo "Installing packages from requirements file..."
    pip3 install -r requirement/requirements_go_analysis.txt
    
    if [ $? -eq 0 ]; then
        echo "✓ All packages installed successfully!"
    else
        echo "✗ Some packages failed to install"
        echo "Trying individual package installation..."
        
        # Install packages individually
        pip3 install pandas numpy matplotlib seaborn scipy statsmodels openpyxl
        
        if [ $? -eq 0 ]; then
            echo "✓ Individual package installation successful!"
        else
            echo "✗ Package installation failed"
            exit 1
        fi
    fi
else
    echo "Requirements file not found, installing packages individually..."
    pip3 install pandas numpy matplotlib seaborn scipy statsmodels openpyxl
    
    if [ $? -eq 0 ]; then
        echo "✓ Individual package installation successful!"
    else
        echo "✗ Package installation failed"
        exit 1
    fi
fi

echo ""
echo "Verifying installation..."
echo "========================"

# Test if packages can be imported
python3 -c "
import sys
packages = ['pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels', 'openpyxl']
failed = []

for package in packages:
    try:
        __import__(package)
        print(f'✓ {package}')
    except ImportError:
        print(f'✗ {package}')
        failed.append(package)

if failed:
    print(f'\nFailed packages: {', '.join(failed)}')
    sys.exit(1)
else:
    print('\n✓ All packages imported successfully!')
"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Dependency installation completed successfully!"
    echo ""
    echo "Next steps:"
    echo "1. Run tests: python tests/test_go_analysis.py"
    echo "2. Run analysis: python src/go_term_analysis.py"
    echo "3. Or use runner script: ./scripts/run_go_analysis.sh"
    echo ""
    
    if [ -d "go_analysis_env" ]; then
        echo "Note: Virtual environment is active. To deactivate, run:"
        echo "  deactivate"
        echo ""
        echo "To reactivate later, run:"
        echo "  source go_analysis_env/bin/activate"
    fi
else
    echo ""
    echo "✗ Some packages failed to import"
    echo "Please check the error messages above and try reinstalling"
    exit 1
fi

echo "Installation script completed!"
