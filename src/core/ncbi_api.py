"""
NCBI API Client - Classe para interação com banco de dados NCBI
"""

from Bio import Entrez
import requests
from typing import List, Dict, Optional


class NCBIClient:
    """Cliente para interação com API do NCBI"""
    
    def __init__(self, email: str):
        """
        Inicializa cliente NCBI
        
        Args:
            email: Email para identificação na API NCBI
        """
        self.email = email
        Entrez.email = email
        
    def search_organism(self, organism_name: str, database: str = 'nucleotide', 
                       max_results: int = 100) -> List[str]:
        """
        Busca organismo no banco de dados NCBI
        
        Args:
            organism_name: Nome do organismo
            database: Banco de dados (nucleotide, genome, taxonomy, etc)
            max_results: Número máximo de resultados
            
        Returns:
            Lista de IDs encontrados
        """
        try:
            search_handle = Entrez.esearch(
                db=database,
                term=organism_name,
                retmax=max_results
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            return search_results.get('IdList', [])
        except Exception as e:
            raise Exception(f"Erro ao buscar {organism_name}: {str(e)}")
    
    def fetch_sequence(self, seq_id: str, database: str = 'nucleotide') -> Dict:
        """
        Busca sequência por ID
        
        Args:
            seq_id: ID da sequência
            database: Banco de dados
            
        Returns:
            Dicionário com dados da sequência
        """
        try:
            fetch_handle = Entrez.efetch(
                db=database,
                id=seq_id,
                rettype='gb',
                retmode='text'
            )
            record = fetch_handle.read()
            fetch_handle.close()
            
            return {
                'id': seq_id,
                'record': record,
                'database': database
            }
        except Exception as e:
            raise Exception(f"Erro ao buscar sequência {seq_id}: {str(e)}")
    
    def fetch_taxonomy(self, organism_name: str) -> Optional[Dict]:
        """
        Busca informações taxonômicas
        
        Args:
            organism_name: Nome do organismo
            
        Returns:
            Dicionário com informações taxonômicas ou None
        """
        try:
            search_handle = Entrez.esearch(db='taxonomy', term=organism_name)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results.get('IdList', [])
            if not id_list:
                return None
            
            tax_id = id_list[0]
            fetch_handle = Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml')
            records = Entrez.read(fetch_handle)
            fetch_handle.close()
            
            if records:
                return records[0]
            return None
        except Exception as e:
            raise Exception(f"Erro ao buscar taxonomia: {str(e)}")
    
    def fetch_genome_info(self, organism_name: str) -> List[Dict]:
        """
        Busca informações de genoma
        
        Args:
            organism_name: Nome do organismo
            
        Returns:
            Lista de dicionários com informações de genoma
        """
        try:
            search_handle = Entrez.esearch(
                db='genome',
                term=organism_name,
                retmax=10
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results.get('IdList', [])
            genomes = []
            
            for genome_id in id_list:
                fetch_handle = Entrez.efetch(
                    db='genome',
                    id=genome_id,
                    retmode='xml'
                )
                record = Entrez.read(fetch_handle)
                fetch_handle.close()
                genomes.append(record)
            
            return genomes
        except Exception as e:
            raise Exception(f"Erro ao buscar genoma: {str(e)}")


# Exemplo de uso (comentado para não executar):
"""
if __name__ == "__main__":
    client = NCBIClient(email="seu.email@exemplo.com")
    
    # Buscar organismo
    ids = client.search_organism("Homo sapiens", max_results=5)
    print(f"IDs encontrados: {ids}")
    
    # Buscar sequência
    if ids:
        seq_data = client.fetch_sequence(ids[0])
        print(f"Sequência: {seq_data['id']}")
    
    # Buscar taxonomia
    tax_info = client.fetch_taxonomy("Homo sapiens")
    if tax_info:
        print(f"Nome científico: {tax_info.get('ScientificName')}")
"""
