# 🧬 Consulta NCBI - DNA e Genoma

Programa de consulta ao banco de dados do NCBI (National Center for Biotechnology Information) para buscar informações completas sobre DNA, genoma, sequências genéticas e dados de organismos.

## 🎯 Arquitetura Modular (OOP)

Este projeto foi refatorado de um arquivo monolítico de 2229 linhas para uma **arquitetura modular orientada a objetos**, melhorando manutenibilidade, testabilidade e legibilidade:

```
📦 Projeto
├── 📄 main.py                 # Ponto de entrada
├── 📁 src/
│   ├── 📄 config.py          # Configurações (UI_COLORS, FONTS, ALIGNMENT_CONFIG)
│   ├── 📄 translations.py    # Suporte multilíngue (7 idiomas)
│   ├── 📁 core/
│   │   ├── 📄 ncbi_api.py    # Classe NCBIClient - API NCBI
│   │   └── 📄 alignment.py   # Classe SequenceAligner - Alinhamento DNA/Proteína
│   ├── 📁 ui/               # Componentes de interface (em desenvolvimento)
│   └── 📁 export/           # Gerenciadores de exportação (em desenvolvimento)
└── 📄 ncbi_dna_search.py     # Arquivo original (mantido como referência)
```

**Benefícios da Refatoração:**
- ✅ Cada arquivo < 500 linhas (vs. 2229 original)
- ✅ Classes com responsabilidade única
- ✅ Testabilidade independente
- ✅ Reutilização de componentes
- ✅ Colaboração facilitada

## 📋 Funcionalidades

- ✅ Busca de **sequências de nucleotídeos** (DNA/RNA)
- ✅ Informações **taxonômicas** completas
- ✅ Dados de **genoma e assemblies**
- ✅ Informações sobre **genes**
- ✅ **Publicações científicas** relacionadas (PubMed)
- ✅ **Comparação de sequências DNA/Proteína** com alinhamento estilo VectorBuilder
- ✅ **Alinhamento proteico com BLOSUM62** e teste de 6 frames
- ✅ **Exportação multilíngue** (PDF/FASTA) em 7 idiomas
- ✅ Interface gráfica simples e eficiente
- ✅ Visualização de sequências genéticas formatadas
- ✅ Múltiplas abas organizadas por tipo de informação

## 🚀 Como Usar

### 1. Instalação das Dependências

Abra o terminal no diretório do projeto e execute:

```bash
pip install -r requirements.txt
```

### 2. Configurar Email

No arquivo [ncbi_dna_search.py](ncbi_dna_search.py#L75), altere o email do NCBI:

```python
Entrez.email = "seu_email@exemplo.com"
```

Coloque seu email real (exigência da NCBI API).

### 3. Executar o Programa

```bash
python main.py
```

Ou diretamente:
```bash
python ncbi_dna_search.py
```

### 4. Realizar Buscas

1. Digite o nome do organismo no campo de busca (ex: "Homo sapiens", "Escherichia coli", "Canis lupus")
2. Clique em **🔍 Buscar** ou pressione Enter
3. Explore os resultados nas abas:
   - **📋 Informações Gerais**: Resumo e dados dos genes
   - **🧬 Sequências**: Sequências de DNA completas
   - **🌳 Taxonomia**: Classificação taxonômica
   - **🔬 Genoma**: Informações do genoma e assemblies
   - **📚 Publicações**: Artigos científicos relacionados

### 5. Comparar Sequências

Na aba **Comparação**, você pode:
- **Alinhamento DNA**: Compara sequências nucleotídicas
- **Alinhamento Proteína**: Traduz DNA→Proteína e alinha com BLOSUM62
  - Testa 6 frames (3 forward + 3 reverse complement)
  - Seleciona automaticamente o melhor frame
  - Score idêntico ao VectorBuilder

### 6. Exportar Resultados

Clique em **💾 Exportar** para salvar em:
- **PDF**: Relatório multilíngue com estatísticas completas
- **FASTA**: Formato padrão para sequências biológicas

## 🧩 Arquitetura Modular - Classes Principais

### NCBIClient (`src/core/ncbi_api.py`)

Gerencia todas as interações com a API NCBI:

```python
from src.core.ncbi_api import NCBIClient

client = NCBIClient(email="seu_email@exemplo.com")

# Buscar organismo
ids = client.search_organism("Homo sapiens", database="nucleotide", max_results=10)

# Obter sequência
sequence_data = client.fetch_sequence(ids[0], database="nucleotide")

# Informações taxonômicas
taxonomy = client.fetch_taxonomy("Canis lupus")
```

### SequenceAligner (`src/core/alignment.py`)

Alinha sequências DNA e proteínas:

```python
from src.core.alignment import SequenceAligner

aligner = SequenceAligner(
    match_score=2,
    mismatch_score=-1,
    gap_open=-2.0,
    gap_extend=-2.0
)

# Alinhamento DNA
dna_result = aligner.align_dna(seq1, seq2, alignment_type='global')

# Alinhamento Proteína (testa 6 frames)
protein_result = aligner.align_protein(dna1, dna2, use_best_frame=True)
```

Veja [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) para mais exemplos.

## 📊 Bancos de Dados Disponíveis

O programa consulta automaticamente múltiplos bancos de dados do NCBI:

- **Nucleotide**: Sequências de DNA e RNA
- **Genome**: Genomas completos
- **Gene**: Informações sobre genes específicos
- **Protein**: Sequências de proteínas
- **Taxonomy**: Classificação taxonômica
- **PubMed**: Publicações científicas
- **Assembly**: Assemblies de genomas

## 🔍 Exemplos de Busca

- `Homo sapiens` - Humano
- `Canis lupus familiaris` - Cachorro
- `Escherichia coli` - Bactéria E. coli
- `BRCA1` - Gene específico
- `Tyrannosaurus rex` - Dinossauro
- `SARS-CoV-2` - Vírus COVID-19

## 🌍 Suporte Multilíngue

Interface disponível em 7 idiomas:

| Idioma | Código | Status |
|--------|--------|--------|
| 🇧🇷 Português | `pt` | ✅ Completo |
| 🇺🇸 English | `en` | ✅ Complete |
| 🇪🇸 Español | `es` | ✅ Completo |
| 🇫🇷 Français | `fr` | ✅ Complet |
| 🇩🇪 Deutsch | `de` | ✅ Vollständig |
| 🇨🇳 中文 | `zh` | ✅ 完成 |
| 🇷🇺 Русский | `ru` | ✅ Завершено |

Altere o idioma na interface gráfica ou via [translations.py](src/translations.py).

## ⚙️ Requisitos

- Python 3.7 ou superior
- Conexão com a internet
- Tkinter (geralmente já incluído no Python)
- Bibliotecas: Biopython, ReportLab (veja [requirements.txt](requirements.txt))

## 📝 Notas Importantes

1. **Limite de Requisições**: A NCBI API tem limite de 3 requisições por segundo sem API key
2. **Email Obrigatório**: Sempre configure seu email no código
3. **Dados Públicos**: Todos os dados são públicos e de livre acesso
4. **Sequências Grandes**: Sequências muito grandes são truncadas na visualização (primeiros 5000 bp)
5. **Alinhamento VectorBuilder**: O alinhamento proteico usa os mesmos parâmetros do VectorBuilder (BLOSUM62, gap penalties -2.0)

## 🛠️ Tecnologias Utilizadas

- **Python 3.11+**: Linguagem principal
- **Biopython 1.81+**: API NCBI, alinhamentos (pairwise2, PairwiseAligner)
- **Tkinter**: Interface gráfica multiplataforma
- **ReportLab**: Geração de PDFs
- **NCBI E-utilities**: API pública do NCBI
- **BLOSUM62**: Matriz de substituição para alinhamento proteico

## 📖 Documentação

- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) - API oficial NCBI
- [Biopython Tutorial](https://biopython.org/wiki/Documentation) - Documentação Biopython
- [ARCHITECTURE.md](ARCHITECTURE.md) - Arquitetura do projeto
- [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) - Guia de migração e exemplos

## 🗺️ Roadmap

### ✅ Fase 1 - Infraestrutura (Completo)
- [x] Estrutura modular de diretórios
- [x] Classes NCBIClient e SequenceAligner
- [x] Sistema de configuração e traduções
- [x] Documentação completa

### 🔄 Fase 2 - UI Modular (Próximo)
- [ ] Extrair componentes de UI do monolito
- [ ] Classes para MainWindow, Tabs, SearchBar
- [ ] Padrão Observer para atualizações

### 📦 Fase 3 - Sistema de Export
- [ ] ExportManager com estratégias
- [ ] PDFExporter e FASTAExporter
- [ ] Templates multilíngues

### 🧪 Fase 4 - Testes
- [ ] Unit tests para NCBIClient
- [ ] Unit tests para SequenceAligner
- [ ] Testes de integração
- [ ] CI/CD pipeline

### 🚀 Fase 5 - Finalização
- [ ] Deprecar arquivo monolítico
- [ ] Otimização de performance
- [ ] Documentação final

## 🤝 Contribuições

Contribuições são bem-vindas! Veja [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) para entender a estrutura do projeto.

**Como contribuir:**
1. Fork o projeto
2. Crie uma branch para sua feature (`git checkout -b feature/MinhaFeature`)
3. Commit suas mudanças (`git commit -m 'Adiciona MinhaFeature'`)
4. Push para a branch (`git push origin feature/MinhaFeature`)
5. Abra um Pull Request

## 📄 Licença

Uso educacional e científico. Respeite os termos de uso da NCBI API.

## 👨‍💻 Autor

Projeto educacional desenvolvido para consulta e análise de dados biológicos do NCBI.
