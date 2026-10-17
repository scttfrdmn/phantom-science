#!/bin/bash

echo "=== INSTALLING GENOMICS COMPARISON TOOLS ==="
echo "Setting up ADMIXTURE, PLINK, and fastSTRUCTURE for benchmarking"
echo ""

# Create tools directory
mkdir -p genomics_tools
cd genomics_tools

echo "📥 Installing PLINK 2.0..."
if ! command -v plink2 &> /dev/null; then
    # Download PLINK2 for macOS ARM64
    wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_mac_arm64_latest.zip
    if [ -f "plink2_mac_arm64_latest.zip" ]; then
        unzip -q plink2_mac_arm64_latest.zip
        chmod +x plink2
        sudo mv plink2 /usr/local/bin/ 2>/dev/null || cp plink2 ~/bin/ 2>/dev/null || echo "PLINK2 downloaded to $(pwd)/plink2"
        echo "✅ PLINK2 installed"
    else
        echo "❌ PLINK2 download failed"
    fi
else
    echo "✅ PLINK2 already installed"
fi

echo ""
echo "📥 Installing ADMIXTURE..."
if ! command -v admixture &> /dev/null; then
    # Download ADMIXTURE for macOS
    wget -q http://dalexander.github.io/admixture/binaries/admixture_macosx-1.3.0.tar.gz
    if [ -f "admixture_macosx-1.3.0.tar.gz" ]; then
        tar -xzf admixture_macosx-1.3.0.tar.gz
        chmod +x admixture_macosx-1.3.0/admixture
        sudo mv admixture_macosx-1.3.0/admixture /usr/local/bin/ 2>/dev/null || cp admixture_macosx-1.3.0/admixture ~/bin/ 2>/dev/null || echo "ADMIXTURE downloaded to $(pwd)/admixture_macosx-1.3.0/admixture"
        echo "✅ ADMIXTURE installed"
    else
        echo "❌ ADMIXTURE download failed"
    fi
else
    echo "✅ ADMIXTURE already installed"  
fi

echo ""
echo "📥 Installing fastSTRUCTURE..."
if ! command -v structure.py &> /dev/null; then
    # Install fastSTRUCTURE via conda/pip if available
    if command -v conda &> /dev/null; then
        conda install -c bioconda faststructure -y
        echo "✅ fastSTRUCTURE installed via conda"
    elif command -v pip3 &> /dev/null; then
        pip3 install numpy scipy cython
        git clone https://github.com/rajanil/fastStructure.git
        cd fastStructure
        python3 setup.py build_ext --inplace
        echo "✅ fastSTRUCTURE built from source"
        cd ..
    else
        echo "⚠️  fastSTRUCTURE requires conda or pip3 - skipping for now"
    fi
else
    echo "✅ fastSTRUCTURE already available"
fi

echo ""
echo "🔧 Setting up PATH..."
export PATH="$(pwd):$(pwd)/admixture_macosx-1.3.0:$(pwd)/fastStructure:$PATH"
echo "export PATH=\"$(pwd):$(pwd)/admixture_macosx-1.3.0:$(pwd)/fastStructure:\$PATH\"" >> ~/.bashrc 2>/dev/null || true
echo "export PATH=\"$(pwd):$(pwd)/admixture_macosx-1.3.0:$(pwd)/fastStructure:\$PATH\"" >> ~/.zshrc 2>/dev/null || true

echo ""
echo "✅ Genomics tools installation completed!"
echo "📁 Tools installed in: $(pwd)"
echo ""
echo "🧬 Available tools:"
command -v plink2 && echo "  ✅ PLINK2: $(command -v plink2)" || echo "  📍 PLINK2: $(pwd)/plink2"
command -v admixture && echo "  ✅ ADMIXTURE: $(command -v admixture)" || echo "  📍 ADMIXTURE: $(pwd)/admixture_macosx-1.3.0/admixture" 
ls fastStructure/structure.py 2>/dev/null && echo "  ✅ fastSTRUCTURE: $(pwd)/fastStructure/structure.py" || echo "  ⚠️  fastSTRUCTURE: Not installed"

echo ""
echo "🎯 Ready for genomic sequence benchmarking!"