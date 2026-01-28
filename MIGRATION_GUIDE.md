# 📊 Guia de Migração - Estrutura Modular

## ✅ O que foi criado

### 1. Estrutura de Diretórios
```
src/
├── config.py              ✅ Criado - Configurações globais
├── translations.py        ✅ Criado - 7 idiomas completos
├── core/
│   ├── __init__.py       ✅ Criado
│   ├── ncbi_api.py       ✅ Criado - Classe NCBIClient
│   └── alignment.py      ✅ Criado - Classe SequenceAligner
├── ui/
│   └── __init__.py       ✅ Criado
└── export/
    └── __init__.py       ✅ Criado
```

### 2. Arquivos na Raiz
```
main.py                   ✅ Criado - Novo ponto de entrada
ARCHITECTURE.md           ✅ Criado - Documentação arquitetura
MIGRATION_GUIDE.md        ✅ Criado - Este guia
ncbi_dna_search.py        ⚠️  Original mantido (funcional)
```

## 🎯 Como Usar Agora

### Opção 1: Arquivo Original (Ainda Funcional)
```bash
python ncbi_dna_search.py
```
✅ Continua funcionando 100%  
✅ Todas as 2229 linhas intactas  
✅ Sem mudanças no comportamento

### Opção 2: Nova Estrutura Modular
```bash
python main.py
```
⚠️  Por enquanto chama o original  
🔄 Migração gradual em andamento

## 📦 Módulos Criados e Prontos para Uso

### 1. **NCBIClient** (src/core/ncbi_api.py)

```python
from src.core.ncbi_api import NCBIClient

# Criar cliente
client = NCBIClient(email="seu@email.com")

# Buscar organismo
ids = client.search_organism("Homo sapiens", max_results=10)

# Buscar sequência
seq_data = client.fetch_sequence(ids[0])

# Buscar taxonomia
tax_info = client.fetch_taxonomy("Homo sapiens")

# Buscar genoma
genomes = client.fetch_genome_info("Homo sapiens")
```

**Métodos disponíveis:**
- `search_organism()` - Busca IDs no NCBI
- `fetch_sequence()` - Busca sequência por ID
- `fetch_taxonomy()` - Informações taxonômicas
- `fetch_genome_info()` - Informações de genoma

### 2. **SequenceAligner** (src/core/alignment.py)

```python
from src.core.alignment import SequenceAligner

# Criar alinhador
aligner = SequenceAligner(
    match_score=2.0,
    mismatch_score=-1.0,
    gap_open=-2.0,
    gap_extend=-2.0
)

# Alinhamento DNA
result = aligner.align_dna(
    seq1="ATCGATCG",
    seq2="ATGGATCG",
    alignment_type='global'  # ou 'local'
)
print(f"Identity: {result['identity']:.2f}%")
print(f"Score: {result['score']}")

# Alinhamento Proteína (com best frame)
result = aligner.align_protein(
    dna1="ATGCGATCGATCG",
    dna2="ATGCGATGGATCG",
    use_best_frame=True  # Testa todos os 6 frames
)
print(f"Frame 1: {result['frame1']}, Frame 2: {result['frame2']}")
print(f"Score: {result['score']}")
```

**Recursos:**
- Alinhamento DNA: Global e Local
- Alinhamento Proteína: Com seleção automática do melhor frame
- Estatísticas detalhadas: Identity, Similarity, Matches, Gaps
- Suporte BLOSUM62 para proteínas

### 3. **Translations** (src/translations.py)

```python
from src.translations import get_translation, TRANSLATIONS

# Usar função helper
texto = get_translation('pt', 'search_button')  # "🔍  Buscar"
texto = get_translation('en', 'search_button')  # "🔍  Search"

# Acessar diretório direto
titulo_pt = TRANSLATIONS['pt']['title']
titulo_es = TRANSLATIONS['es']['title']
```

**Idiomas disponíveis:**
- 🇧🇷 Português (pt)
- 🇺🇸 English (en)
- 🇪🇸 Español (es)
- 🇫🇷 Français (fr)
- 🇩🇪 Deutsch (de)
- 🇨🇳 中文 (zh)
- 🇷🇺 Русский (ru)

### 4. **Config** (src/config.py)

```python
from src.config import UI_COLORS, FONTS, ALIGNMENT_CONFIG

# Cores padronizadas
window.configure(bg=UI_COLORS['bg_dark'])
button.configure(bg=UI_COLORS['accent_red'])

# Fontes padronizadas
label.configure(font=FONTS['heading'])

# Configurações de alinhamento
match_score = ALIGNMENT_CONFIG['match_score']
```

## 🔄 Próximos Passos da Migração

### Fase 2: UI Modules (Próxima)
- [ ] Criar `src/ui/main_window.py`
- [ ] Extrair componentes de UI
- [ ] Separar tabs em arquivos individuais
- [ ] Implementar padrão Observer para atualização

### Fase 3: Export Module
- [ ] Criar `src/export/export_manager.py`
- [ ] Classe `PDFExporter`
- [ ] Classe `FASTAExporter`
- [ ] Suporte multilíngue automático

### Fase 4: Testes
- [ ] Criar `tests/` directory
- [ ] Testes unitários para NCBIClient
- [ ] Testes unitários para SequenceAligner
- [ ] Testes de integração

### Fase 5: Documentação
- [ ] Docstrings completos
- [ ] Exemplos de uso
- [ ] Tutorial passo a passo
- [ ] API Reference

## 💡 Benefícios Imediatos

### Para Desenvolvimento
- ✅ Código organizado por responsabilidade
- ✅ Fácil de entender cada módulo
- ✅ Reutilização de código
- ✅ Testes independentes possíveis

### Para Manutenção
- ✅ Mudanças localizadas
- ✅ Menos risco de quebrar outras partes
- ✅ Mais fácil de debugar
- ✅ Melhor rastreabilidade

### Para Colaboração
- ✅ Múltiplas pessoas podem trabalhar simultaneamente
- ✅ Pull requests menores e focados
- ✅ Code review mais eficiente
- ✅ Onboarding mais rápido

## 🚀 Como Contribuir Agora

### 1. Escolha um Módulo
Escolha uma parte para refatorar:
- UI Components
- Export Manager
- Utils & Helpers
- Documentation

### 2. Siga o Padrão
```python
"""
Docstring descritivo do módulo
"""

class MyClass:
    """Docstring da classe"""
    
    def __init__(self, param: type):
        """
        Inicializador
        
        Args:
            param: Descrição do parâmetro
        """
        self.param = param
    
    def method(self) -> return_type:
        """
        Descrição do método
        
        Returns:
            Descrição do retorno
        """
        pass
```

### 3. Teste Individualmente
```python
if __name__ == "__main__":
    # Código de teste aqui
    obj = MyClass(param="test")
    result = obj.method()
    print(f"Result: {result}")
```

### 4. Documente
- Adicione docstrings
- Atualize ARCHITECTURE.md
- Adicione exemplos de uso

## 📞 Contato

Para dúvidas sobre a migração, abra uma issue ou contacte o time de desenvolvimento.

---

**Status**: ✅ Fase 1 Completa - Estrutura base criada  
**Próximo**: 🔄 Fase 2 - UI Modules  
**Data**: Janeiro 2026
