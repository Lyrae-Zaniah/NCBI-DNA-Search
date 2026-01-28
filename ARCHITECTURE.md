# NCBI DNA Search - Arquitetura Modular

## 📁 Estrutura do Projeto Refatorado

```
Consulta Em Banco de Dados de DNA/
├── main.py                    # Ponto de entrada principal
├── ncbi_dna_search.py        # Arquivo original (será migrado gradualmente)
├── requirements.txt
├── README.md
└── src/                      # Código fonte modularizado
    ├── config.py             # ⚙️  Configurações globais
    ├── translations.py       # 🌍 Traduções multilíngue
    ├── core/                 # 🧬 Lógica de negócio
    │   ├── __init__.py
    │   ├── ncbi_api.py      # API do NCBI
    │   └── alignment.py     # Alinhamento de sequências
    ├── ui/                   # 🖥️  Interface gráfica
    │   ├── __init__.py
    │   └── main_window.py   # Janela principal
    └── export/               # 💾 Exportação
        ├── __init__.py
        └── export_manager.py # Gerenciamento de exports
```

## 🎯 Benefícios da Arquitetura Modular

### ✅ Antes (Arquivo Único)
- ❌ Mais de 2200 linhas em um arquivo
- ❌ Difícil manutenção
- ❌ Testes complicados
- ❌ Reutilização limitada

### ✅ Depois (Modular)
- ✅ Arquivos pequenos e focados (<500 linhas cada)
- ✅ Separação clara de responsabilidades
- ✅ Fácil de testar individualmente
- ✅ Código reutilizável
- ✅ Melhor organização em equipe

## 📦 Módulos Criados

### 1. `src/config.py`
- Configurações globais da aplicação
- Cores da UI
- Parâmetros de alinhamento
- Configurações SSL/NCBI

### 2. `src/translations.py`
- Dicionários de tradução para 7 idiomas
- Função helper `get_translation()`
- Português, English, Español, Français, Deutsch, 中文, Русский

### 3. `src/core/` (A ser implementado)
- **ncbi_api.py**: Integração com API NCBI
  - Classe `NCBIClient`
  - Métodos para buscar sequências, taxonomia, genoma
  
- **alignment.py**: Algoritmos de alinhamento
  - Classe `SequenceAligner`
  - Alinhamento DNA/Proteína
  - Cálculo de similaridade

### 4. `src/ui/` (A ser implementado)
- **main_window.py**: Interface gráfica principal
  - Classe `MainWindow`
  - Gerenciamento de tabs
  - Componentes reutilizáveis

### 5. `src/export/` (A ser implementado)
- **export_manager.py**: Exportação de dados
  - Classe `ExportManager`
  - Export para PDF, FASTA
  - Formatação multilíngue

## 🚀 Como Executar

```bash
# Método 1: Usando o novo main.py
python main.py

# Método 2: Arquivo original (ainda funcional)
python ncbi_dna_search.py
```

## 📝 Roadmap de Migração

### Fase 1: Estrutura Base ✅ (COMPLETO)
- [x] Criar diretórios src/
- [x] Extrair config.py
- [x] Extrair translations.py
- [x] Criar main.py

### Fase 2: Core Modules (PRÓXIMO)
- [ ] Implementar ncbi_api.py
- [ ] Implementar alignment.py
- [ ] Criar testes unitários

### Fase 3: UI Modules
- [ ] Refatorar main_window.py
- [ ] Separar componentes de UI
- [ ] Implementar padrão Observer

### Fase 4: Export & Utils
- [ ] Implementar export_manager.py
- [ ] Adicionar logging
- [ ] Documentação completa

### Fase 5: Finalização
- [ ] Migração completa
- [ ] Deprecar ncbi_dna_search.py
- [ ] Testes de integração

## 🔧 Padrões de Design Utilizados

- **MVC (Model-View-Controller)**
  - Model: `src/core/`
  - View: `src/ui/`
  - Controller: Conexão entre ambos

- **Singleton**: Para configurações globais
- **Factory**: Para criação de exportadores
- **Strategy**: Para diferentes tipos de alinhamento

## 📚 Dependências

```
biopython
tkinter
matplotlib
reportlab
certifi
```

## 👥 Contribuindo

Com a estrutura modular, contribuir ficou mais fácil:

1. Escolha um módulo específico
2. Faça alterações em arquivos pequenos
3. Teste individualmente
4. Pull request focado

## 📄 Licença

[Sua licença aqui]

---

**Nota**: A migração está em andamento. O arquivo original `ncbi_dna_search.py` ainda é funcional e pode ser usado normalmente. A nova estrutura modular será adotada gradualmente.
