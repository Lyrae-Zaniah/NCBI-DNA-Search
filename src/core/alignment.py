"""
Sequence Alignment - Classes para alinhamento de sequências DNA/Proteína
"""

from Bio import Align
from Bio.Seq import Seq
from Bio import pairwise2
from typing import Tuple, Dict, Optional
from collections import Counter


class SequenceAligner:
    """Classe para alinhamento de sequências"""
    
    def __init__(self, match_score: float = 2.0, mismatch_score: float = -1.0,
                 gap_open: float = -2.0, gap_extend: float = -2.0):
        """
        Inicializa alinhador com parâmetros
        
        Args:
            match_score: Score para match
            mismatch_score: Score para mismatch  
            gap_open: Penalidade para abrir gap
            gap_extend: Penalidade para extender gap
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
    def align_dna(self, seq1: str, seq2: str, 
                  alignment_type: str = 'global') -> Dict:
        """
        Alinha duas sequências de DNA
        
        Args:
            seq1: Primeira sequência
            seq2: Segunda sequência
            alignment_type: 'global' ou 'local'
            
        Returns:
            Dicionário com resultados do alinhamento
        """
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        
        if alignment_type == 'local':
            alignments = pairwise2.align.localms(
                seq1, seq2,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend
            )
        else:
            alignments = pairwise2.align.globalms(
                seq1, seq2,
                self.match_score, self.mismatch_score,
                self.gap_open, self.gap_extend
            )
        
        if not alignments:
            return {'error': 'Nenhum alinhamento encontrado'}
        
        best_alignment = alignments[0]
        aligned_seq1, aligned_seq2, score, begin, end = best_alignment
        
        # Calcula estatísticas
        stats = self._calculate_dna_statistics(aligned_seq1, aligned_seq2)
        
        return {
            'seq1': aligned_seq1,
            'seq2': aligned_seq2,
            'score': score,
            'begin': begin,
            'end': end,
            'type': alignment_type,
            **stats
        }
    
    def align_protein(self, dna1: str, dna2: str, 
                     use_best_frame: bool = True) -> Dict:
        """
        Alinha sequências traduzindo DNA para proteína
        
        Args:
            dna1: Primeira sequência DNA
            dna2: Segunda sequência DNA
            use_best_frame: Se True, testa todos os frames e escolhe melhor
            
        Returns:
            Dicionário com resultados do alinhamento
        """
        if use_best_frame:
            return self._align_protein_best_frame(dna1, dna2)
        else:
            return self._align_protein_single_frame(dna1, dna2, frame=0)
    
    def _align_protein_best_frame(self, dna1: str, dna2: str) -> Dict:
        """Testa todos os frames e retorna melhor alinhamento"""
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        
        best_score = float('-inf')
        best_result = None
        
        # Testa 6 combinações de frames
        for frame1 in range(3):
            seq1_frame = dna1[frame1:]
            seq1_frame = seq1_frame[:len(seq1_frame) - (len(seq1_frame) % 3)]
            protein1 = str(Seq(seq1_frame).translate(to_stop=False))
            
            for frame2 in range(3):
                seq2_frame = dna2[frame2:]
                seq2_frame = seq2_frame[:len(seq2_frame) - (len(seq2_frame) % 3)]
                protein2 = str(Seq(seq2_frame).translate(to_stop=False))
                
                try:
                    alignment = next(aligner.align(protein1, protein2))
                    if alignment.score > best_score:
                        best_score = alignment.score
                        stats = self._calculate_protein_statistics(
                            str(alignment[0]), str(alignment[1])
                        )
                        best_result = {
                            'seq1': str(alignment[0]),
                            'seq2': str(alignment[1]),
                            'score': alignment.score,
                            'frame1': frame1 + 1,
                            'frame2': frame2 + 1,
                            'protein1_length': len(protein1),
                            'protein2_length': len(protein2),
                            **stats
                        }
                except StopIteration:
                    continue
        
        return best_result if best_result else {'error': 'Nenhum alinhamento encontrado'}
    
    def _align_protein_single_frame(self, dna1: str, dna2: str, 
                                    frame: int = 0) -> Dict:
        """Alinha proteínas usando frame específico"""
        # Traduz DNA para proteína
        protein1 = str(Seq(dna1).translate(to_stop=False))
        protein2 = str(Seq(dna2).translate(to_stop=False))
        
        # Configura alinhador
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        
        alignment = next(aligner.align(protein1, protein2))
        
        stats = self._calculate_protein_statistics(
            str(alignment[0]), str(alignment[1])
        )
        
        return {
            'seq1': str(alignment[0]),
            'seq2': str(alignment[1]),
            'score': alignment.score,
            'frame': frame + 1,
            **stats
        }
    
    def _calculate_dna_statistics(self, seq1: str, seq2: str) -> Dict:
        """Calcula estatísticas para alinhamento DNA"""
        matches = 0
        transitions = 0  # A↔G, C↔T
        transversions = 0
        gaps = 0
        
        for a, b in zip(seq1, seq2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                purines = ['A', 'G']
                pyrimidines = ['C', 'T']
                if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                    transitions += 1
                else:
                    transversions += 1
        
        length = len(seq1)
        identity = (matches / length * 100) if length > 0 else 0
        similarity = ((matches + transitions) / length * 100) if length > 0 else 0
        
        return {
            'identity': identity,
            'similarity': similarity,
            'matches': matches,
            'transitions': transitions,
            'transversions': transversions,
            'gaps': gaps,
            'length': length
        }
    
    def _calculate_protein_statistics(self, seq1: str, seq2: str) -> Dict:
        """Calcula estatísticas para alinhamento proteína"""
        matches = 0
        similar = 0
        gaps = 0
        
        # Grupos de aminoácidos similares
        similar_groups = [
            set('RKH'),    # Básicos
            set('DE'),     # Ácidos
            set('AVLIM'),  # Hidrofóbicos
            set('FYW'),    # Aromáticos
            set('STNQ')    # Polares
        ]
        
        for a, b in zip(seq1, seq2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                # Verifica se são similares
                for group in similar_groups:
                    if a in group and b in group:
                        similar += 1
                        break
        
        length = len(seq1)
        identity = (matches / length * 100) if length > 0 else 0
        similarity = ((matches + similar) / length * 100) if length > 0 else 0
        
        return {
            'identity': identity,
            'similarity': similarity,
            'matches': matches,
            'similar': similar,
            'gaps': gaps,
            'length': length
        }


# Exemplo de uso (comentado):
"""
if __name__ == "__main__":
    aligner = SequenceAligner()
    
    # Alinhamento DNA
    seq1 = "ATCGATCG"
    seq2 = "ATGGATCG"
    result = aligner.align_dna(seq1, seq2, 'global')
    print(f"Score: {result['score']}, Identity: {result['identity']:.2f}%")
    
    # Alinhamento proteína
    dna1 = "ATGCGATCGATCG"
    dna2 = "ATGCGATGGATCG"
    result = aligner.align_protein(dna1, dna2)
    print(f"Score: {result['score']}, Frames: {result['frame1']}, {result['frame2']}")
"""
