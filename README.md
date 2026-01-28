# 🧬 NCBI DNA Search

Desktop application for querying NCBI (National Center for Biotechnology Information) databases to retrieve DNA sequences, genome data, and organism information.

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![BioPython](https://img.shields.io/badge/BioPython-1.81+-orange.svg)](https://biopython.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## ✨ Features

- 🔍 **Search** nucleotide sequences (DNA/RNA), genes, genomes, and publications
- 🧬 **Sequence alignment** - DNA and Protein (BLOSUM62 matrix, 6-frame testing)
- 📊 **Statistical analysis** - GC/AT content, base composition, interactive charts
- 🌍 **Multilingual** - 7 languages (EN, PT, ES, FR, DE, ZH, RU)
- 💾 **Export** - PDF reports and FASTA format
- 🎨 **Modern GUI** - Dark theme with organized tabs

## 🚀 Quick Start

```bash
# Clone repository
git clone https://github.com/your-username/ncbi-dna-search.git
cd ncbi-dna-search

# Install dependencies
pip install -r requirements.txt

# Run application
python ncbi_dna_search.py
```

**Note**: Configure your email in the app (NCBI requirement)

## � Alignment Capabilities

### DNA Alignment
- **Global**: End-to-end alignment (same species)
- **Local**: Find regions of similarity (different species)

### Protein Alignment (VectorBuilder-style)
- Translates DNA → Protein automatically
- Tests all 6 reading frames (3 forward + 3 reverse)
- Uses BLOSUM62 substitution matrix
- Calculates: Identity, Similarity, Gaps, Transitions/Transversions

**Example Results:**
```
DNA: 5040 bp → Protein: 1680 aa
Score: 509.0
Identity: 17.03% | Similarity: 24.56% | Gaps: 57.24%
```

## 🧩 Modular Architecture

Refactored from monolithic 2229-line file to modular OOP design:

```
src/
├── config.py          # UI colors, fonts, alignment settings
├── translations.py    # 7 language translations
└── core/
    ├── ncbi_api.py    # NCBIClient class
    └── alignment.py   # SequenceAligner class
```

**Usage Example:**
```python
from src.core.ncbi_api import NCBIClient
from src.core.alignment import SequenceAligner

# NCBI queries
client = NCBIClient(email="your@email.com")
ids = client.search_organism("Homo sapiens")
seq = client.fetch_sequence(ids[0])

# Sequence alignment
aligner = SequenceAligner()
result = aligner.align_protein(dna1, dna2, use_best_frame=True)
```

See [ARCHITECTURE.md](ARCHITECTURE.md) and [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for details.

## 🌍 Languages

| Language | Code | | Language | Code |
|----------|------|---|----------|------|
| 🇺🇸 English | `en` | | 🇫🇷 Français | `fr` |
| 🇧🇷 Português | `pt` | | 🇩🇪 Deutsch | `de` |
| 🇪🇸 Español | `es` | | 🇨🇳 中文 | `zh` |
| | | | 🇷🇺 Русский | `ru` |

## 📦 Dependencies

```
biopython>=1.81
matplotlib>=3.5.0
reportlab>=3.6.0
certifi>=2021.10.8
```

## 🛠️ Tech Stack

- **Python 3.8+** - Core language
- **BioPython** - NCBI API, alignments (pairwise2, PairwiseAligner)
- **Tkinter** - Cross-platform GUI
- **Matplotlib** - Interactive charts
- **ReportLab** - PDF generation
- **BLOSUM62** - Protein substitution matrix

## 📖 Documentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - Project architecture and roadmap
- [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) - Migration guide with code examples
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) - Official NCBI API docs

## 🤝 Contributing

Contributions welcome! See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for project structure.

1. Fork the repository
2. Create feature branch (`git checkout -b feature/NewFeature`)
3. Commit changes (`git commit -m 'Add NewFeature'`)
4. Push to branch (`git push origin feature/NewFeature`)
5. Open Pull Request

## 📄 License

Educational and scientific use. Please respect NCBI API terms of service.

---

**Made with ❤️ and 🧬** | [Report Issues](../../issues) | [Documentation](ARCHITECTURE.md)
