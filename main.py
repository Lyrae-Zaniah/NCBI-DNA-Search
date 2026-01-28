"""
NCBI DNA Search Application - Refatorado com Arquitetura Modular

Estrutura do Projeto:
- src/
  ├── config.py          - Configurações globais
  ├── translations.py    - Dicionários de tradução multilíngue
  ├── core/              - Lógica de negócio
  │   ├── ncbi_api.py    - Integração com NCBI API
  │   └── alignment.py   - Algoritmos de alinhamento
  ├── ui/                - Interface gráfica
  │   └── main_window.py - Janela principal
  └── export/            - Exportação de dados
      └── export_manager.py - Gerenciamento de exportações

Executar: python main.py
"""

import sys
import os

# Adiciona o diretório src ao path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Importa o arquivo original (enquanto não terminamos a refatoração completa)
# TODO: Substituir pela versão modular
import ncbi_dna_search

if __name__ == "__main__":
    ncbi_dna_search.main()
