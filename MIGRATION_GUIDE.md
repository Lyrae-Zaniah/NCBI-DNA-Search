# 📊 Migration Guide - Modular Structure

## ✅ What Was Created

### 1. Directory Structure
```
src/
├── config.py              ✅ Created - Global configurations
├── translations.py        ✅ Created - 7 complete languages
├── core/
│   ├── __init__.py       ✅ Created
│   ├── ncbi_api.py       ✅ Created - NCBIClient class
│   └── alignment.py      ✅ Created - SequenceAligner class
├── ui/
│   └── __init__.py       ✅ Created
└── export/
    └── __init__.py       ✅ Created
```

### 2. Root Files
```
main.py                   ✅ Created - New entry point
ARCHITECTURE.md           ✅ Created - Architecture documentation
MIGRATION_GUIDE.md        ✅ Created - This guide
ncbi_dna_search.py        ⚠️  Original kept (functional)
```

## 🎯 How to Use Now

### Option 1: Original File (Still Functional)
```bash
python ncbi_dna_search.py
```
✅ Still works 100%  
✅ All 2229 lines intact  
✅ No behavior changes

### Option 2: New Modular Structure
```bash
python main.py
```
⚠️  Currently calls original  
🔄 Gradual migration in progress

## 📦 Created Modules Ready to Use

### 1. **NCBIClient** (src/core/ncbi_api.py)

```python
from src.core.ncbi_api import NCBIClient

# Create client
client = NCBIClient(email="your@email.com")

# Search organism
ids = client.search_organism("Homo sapiens", max_results=10)

# Fetch sequence
seq_data = client.fetch_sequence(ids[0])

# Fetch taxonomy
tax_info = client.fetch_taxonomy("Homo sapiens")

# Fetch genome
genomes = client.fetch_genome_info("Homo sapiens")
```

**Available Methods:**
- `search_organism()` - Search IDs in NCBI
- `fetch_sequence()` - Fetch sequence by ID
- `fetch_taxonomy()` - Taxonomic information
- `fetch_genome_info()` - Genome information

### 2. **SequenceAligner** (src/core/alignment.py)

```python
from src.core.alignment import SequenceAligner

# Create aligner
aligner = SequenceAligner(
    match_score=2.0,
    mismatch_score=-1.0,
    gap_open=-2.0,
    gap_extend=-2.0
)

# DNA Alignment
result = aligner.align_dna(
    seq1="ATCGATCG",
    seq2="ATGGATCG",
    alignment_type='global'  # or 'local'
)
print(f"Identity: {result['identity']:.2f}%")
print(f"Score: {result['score']}")

# Protein Alignment (with best frame)
result = aligner.align_protein(
    dna1="ATGCGATCGATCG",
    dna2="ATGCGATGGATCG",
    use_best_frame=True  # Tests all 6 frames
)
print(f"Frame 1: {result['frame1']}, Frame 2: {result['frame2']}")
print(f"Score: {result['score']}")
```

**Features:**
- DNA Alignment: Global and Local
- Protein Alignment: Automatic best frame selection
- Detailed Statistics: Identity, Similarity, Matches, Gaps
- BLOSUM62 support for proteins

### 3. **Translations** (src/translations.py)

```python
from src.translations import get_translation, TRANSLATIONS

# Use helper function
text = get_translation('pt', 'search_button')  # "🔍  Buscar"
text = get_translation('en', 'search_button')  # "🔍  Search"

# Direct dictionary access
title_pt = TRANSLATIONS['pt']['title']
title_es = TRANSLATIONS['es']['title']
```

**Available Languages:**
- 🇧🇷 Portuguese (pt)
- 🇺🇸 English (en)
- 🇪🇸 Spanish (es)
- 🇫🇷 French (fr)
- 🇩🇪 German (de)
- 🇨🇳 Chinese (zh)
- 🇷🇺 Russian (ru)

### 4. **Config** (src/config.py)

```python
from src.config import UI_COLORS, FONTS, ALIGNMENT_CONFIG

# Standardized colors
window.configure(bg=UI_COLORS['bg_dark'])
button.configure(bg=UI_COLORS['accent_red'])

# Standardized fonts
label.configure(font=FONTS['heading'])

# Alignment configurations
match_score = ALIGNMENT_CONFIG['match_score']
```

## 🔄 Next Migration Steps

### Phase 2: UI Modules (Next)
- [ ] Create `src/ui/main_window.py`
- [ ] Extract UI components
- [ ] Separate tabs into individual files
- [ ] Implement Observer pattern for updates

### Phase 3: Export Module
- [ ] Create `src/export/export_manager.py`
- [ ] `PDFExporter` class
- [ ] `FASTAExporter` class
- [ ] Automatic multilingual support

### Phase 4: Testing
- [ ] Create `tests/` directory
- [ ] Unit tests for NCBIClient
- [ ] Unit tests for SequenceAligner
- [ ] Integration tests

### Phase 5: Documentation
- [ ] Complete docstrings
- [ ] Usage examples
- [ ] Step-by-step tutorial
- [ ] API Reference

## 💡 Immediate Benefits

### For Development
- ✅ Code organized by responsibility
- ✅ Easy to understand each module
- ✅ Code reusability
- ✅ Independent testing possible

### For Maintenance
- ✅ Localized changes
- ✅ Less risk of breaking other parts
- ✅ Easier to debug
- ✅ Better traceability

### For Collaboration
- ✅ Multiple people can work simultaneously
- ✅ Smaller, focused pull requests
- ✅ More efficient code review
- ✅ Faster onboarding

## 🚀 How to Contribute Now

### 1. Choose a Module
Select a part to refactor:
- UI Components
- Export Manager
- Utils & Helpers
- Documentation

### 2. Follow the Pattern
```python
"""
Descriptive module docstring
"""

class MyClass:
    """Class docstring"""
    
    def __init__(self, param: type):
        """
        Initializer
        
        Args:
            param: Parameter description
        """
        self.param = param
    
    def method(self) -> return_type:
        """
        Method description
        
        Returns:
            Return value description
        """
        pass
```

### 3. Test Individually
```python
if __name__ == "__main__":
    # Test code here
    obj = MyClass(param="test")
    result = obj.method()
    print(f"Result: {result}")
```

### 4. Document
- Add docstrings
- Update ARCHITECTURE.md
- Add usage examples

## 📞 Contact

For questions about the migration, open an issue or contact the development team.

---

**Status**: ✅ Phase 1 Complete - Base structure created  
**Next**: 🔄 Phase 2 - UI Modules  
**Date**: January 2026
