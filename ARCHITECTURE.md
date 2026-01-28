# NCBI DNA Search - Modular Architecture

## 📁 Refactored Project Structure

```
Consulta Em Banco de Dados de DNA/
├── main.py                    # Main entry point
├── ncbi_dna_search.py        # Original file (gradual migration)
├── requirements.txt
├── README.md
└── src/                      # Modularized source code
    ├── config.py             # ⚙️  Global configurations
    ├── translations.py       # 🌍 Multilingual translations
    ├── core/                 # 🧬 Business logic
    │   ├── __init__.py
    │   ├── ncbi_api.py      # NCBI API
    │   └── alignment.py     # Sequence alignment
    ├── ui/                   # 🖥️  Graphical interface
    │   ├── __init__.py
    │   └── main_window.py   # Main window
    └── export/               # 💾 Export functionality
        ├── __init__.py
        └── export_manager.py # Export management
```

## 🎯 Benefits of Modular Architecture

### ✅ Before (Single File)
- ❌ Over 2200 lines in one file
- ❌ Difficult maintenance
- ❌ Complex testing
- ❌ Limited reusability

### ✅ After (Modular)
- ✅ Small, focused files (<500 lines each)
- ✅ Clear separation of concerns
- ✅ Easy to test individually
- ✅ Reusable code
- ✅ Better team organization

## 📦 Created Modules

### 1. `src/config.py`
- Global application settings
- UI colors
- Alignment parameters
- SSL/NCBI configuration

### 2. `src/translations.py`
- Translation dictionaries for 7 languages
- Helper function `get_translation()`
- Portuguese, English, Spanish, French, German, Chinese, Russian

### 3. `src/core/` (To be implemented)
- **ncbi_api.py**: NCBI API integration
  - `NCBIClient` class
  - Methods to fetch sequences, taxonomy, genome
  
- **alignment.py**: Alignment algorithms
  - `SequenceAligner` class
  - DNA/Protein alignment
  - Similarity calculation

### 4. `src/ui/` (To be implemented)
- **main_window.py**: Main graphical interface
  - `MainWindow` class
  - Tab management
  - Reusable components

### 5. `src/export/` (To be implemented)
- **export_manager.py**: Data export
  - `ExportManager` class
  - Export to PDF, FASTA
  - Multilingual formatting

## 🚀 How to Run

```bash
# Method 1: Using new main.py
python main.py

# Method 2: Original file (still functional)
python ncbi_dna_search.py
```

## 📝 Migration Roadmap

### Phase 1: Base Structure ✅ (COMPLETE)
- [x] Create src/ directories
- [x] Extract config.py
- [x] Extract translations.py
- [x] Create main.py

### Phase 2: Core Modules (NEXT)
- [ ] Implement ncbi_api.py
- [ ] Implement alignment.py
- [ ] Create unit tests

### Phase 3: UI Modules
- [ ] Refactor main_window.py
- [ ] Separate UI components
- [ ] Implement Observer pattern

### Phase 4: Export & Utils
- [ ] Implement export_manager.py
- [ ] Add logging
- [ ] Complete documentation

### Phase 5: Finalization
- [ ] Complete migration
- [ ] Deprecate ncbi_dna_search.py
- [ ] Integration tests

## 🔧 Design Patterns Used

- **MVC (Model-View-Controller)**
  - Model: `src/core/`
  - View: `src/ui/`
  - Controller: Connection between both

- **Singleton**: For global configurations
- **Factory**: For creating exporters
- **Strategy**: For different alignment types

## 📚 Dependencies

```
biopython
tkinter
matplotlib
reportlab
certifi
```

## 👥 Contributing

With the modular structure, contributing is easier:

1. Choose a specific module
2. Make changes to small files
3. Test individually
4. Focused pull request

## 📄 License

[Your license here]

---

**Note**: Migration is in progress. The original file `ncbi_dna_search.py` is still functional and can be used normally. The new modular structure will be adopted gradually.
