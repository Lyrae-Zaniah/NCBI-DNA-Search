"""
Programa de Consulta ao Banco de Dados NCBI
Busca informações sobre DNA, genoma, sequências genéticas e dados de organismos
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
from Bio import Entrez, SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import requests
import json
from datetime import datetime
import threading
import ssl
import certifi
from collections import Counter
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
from reportlab.lib import colors
import io

# Email será configurado pelo usuário na interface

# Configuração SSL para evitar erros de certificado
try:
    import urllib.request
    import certifi
    # Cria um contexto SSL com os certificados do certifi
    ssl_context = ssl.create_default_context(cafile=certifi.where())
    # Configura o opener padrão para usar o contexto SSL
    https_handler = urllib.request.HTTPSHandler(context=ssl_context)
    opener = urllib.request.build_opener(https_handler)
    urllib.request.install_opener(opener)
except:
    # Se falhar, desabilita verificação SSL (menos seguro mas funciona)
    ssl._create_default_https_context = ssl._create_unverified_context

# Dicionário de traduções
TRANSLATIONS = {
    'pt': {
        'title': 'Consulta NCBI - DNA e Genoma',
        'email_label': 'Seu Email (obrigatório):',
        'email_placeholder': 'seu.email@exemplo.com',
        'email_required': 'Por favor, insira seu email antes de buscar!',
        'email_invalid': 'Por favor, insira um email válido!',
        'search_label': 'Pesquisar Organismo:',
        'search_button': '🔍  Buscar',
        'clear_button': '🗑  Limpar',
        'language_label': 'Idioma:',
        'ready': 'Pronto para buscar',
        'results_cleared': 'Resultados limpos',
        'searching': '🔍 Buscando informações sobre',
        'completed': '✅ Busca concluída para',
        'error': '❌ Erro na busca:',
        'attention': 'Atenção',
        'enter_organism': 'Digite o nome de um organismo para buscar!',
        'error_title': 'Erro',
        'error_searching': 'Erro ao buscar dados:\n',
        'tab_info': '📋  Informações Gerais',
        'tab_sequences': '🧬  Sequências',
        'tab_taxonomy': '🌳  Taxonomia',
        'tab_genome': '🔬  Genoma',
        'tab_publications': '📚  Publicações',
        'searching_nucleotides': 'BUSCANDO SEQUÊNCIAS DE NUCLEOTÍDEOS PARA:',
        'total_sequences': '📊 Total de sequências encontradas:',
        'showing_first': '📥 Mostrando primeiros',
        'results': 'resultados',
        'register': 'REGISTRO #',
        'annotations': '🔬 Anotações:',
        'sequence': 'SEQUÊNCIA #',
        'size': 'Tamanho:',
        'view_online': '🌐 Ver online:',
        'sequence_unavailable': '⚠️ Sequência não disponível (registro CON/scaffold sem sequência direta)',
        'features': '🧬 Features (primeiros 20):',
        'no_sequences': 'Nenhuma sequência encontrada.',
        'view_browser': '🌐 VISUALIZAR NO NAVEGADOR NCBI',
        'sequences_link': '• Sequências:',
        'genome_complete': '• Genoma Completo:',
        'graphics': '• Graphics:',
        'tax_info': 'INFORMAÇÕES TAXONÔMICAS PARA:',
        'organism': 'ORGANISMO #',
        'scientific_name': 'Nome Científico:',
        'common_name': 'Nome Comum:',
        'other_names': 'Outros Nomes:',
        'genetic_code': '🧬 Código Genético:',
        'no_tax_info': 'Nenhuma informação taxonômica encontrada.',
        'genome_info': 'INFORMAÇÕES DE GENOMA PARA:',
        'total_genomes': '📊 Total de genomas encontrados:',
        'genome': 'GENOMA #',
        'no_genomes': 'Nenhum genoma encontrado.',
        'assemblies': 'ASSEMBLIES DE GENOMA',
        'assembly': 'ASSEMBLY #',
        'view_graphics': '🖼️  Ver Gráficos:',
        'visualization_tools': '🌐 FERRAMENTAS DE VISUALIZAÇÃO NCBI',
        'visual_resources': '📊 RECURSOS VISUAIS:',
        'tip': '💡 Dica: Clique nos links acima para ver mapas, gráficos e imagens interativas!',
        'gene_info': 'INFORMAÇÕES DE GENES',
        'total_genes': '📊 Total de genes encontrados:',
        'gene': 'GENE #',
        'publications': 'PUBLICAÇÕES RELACIONADAS A:',
        'total_publications': '📊 Total de publicações encontradas:',
        'publication': 'PUBLICAÇÃO #',
        'authors': 'Autores:',
        'summary': 'Resumo:',
        'no_publications': 'Nenhuma publicação encontrada.',
        'name': 'Nome:',
        'description': 'Descrição:',
        'division': 'Divisão:',
        'rank': 'Rank:',
        'lineage': '🌳 Linhagem Taxonômica:',
        'organism_label': 'Organismo:',
        'status': 'Status:',
        'access': 'Acesso:',
        'coverage': 'Cobertura:',
        'chromosome': 'Cromossomo:',
        'location': 'Localização:',
        'title_label': 'Título:',
        'year': 'Ano:',
        'gene_id': 'Gene ID:',
        'tab_analysis': '📊  Análise',
        'tab_compare': '🔬  Comparar',
        'export_button': '💾  Exportar',
        'analyze_button': '📊  Analisar',
        'compare_button': '🔍  Comparar',
        'nucleotide_composition': 'COMPOSIÇÃO DE NUCLEOTÍDEOS',
        'gc_content': 'Conteúdo GC:',
        'at_content': 'Conteúdo AT:',
        'total_bases': 'Total de Bases:',
        'sequence_1': 'Sequência 1:',
        'sequence_2': 'Sequência 2:',
        'paste_sequence': 'Cole a sequência aqui...',
        'similarity': 'Similaridade:',
        'alignment_score': 'Score de Alinhamento:',
        'identical_bases': 'Bases Idênticas:',
        'export_success': 'Exportação realizada com sucesso!',
        'success_title': 'Sucesso',
        'error_export': 'Erro ao exportar',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Escolha o formato:',
        'select_sequence': 'Selecione uma sequência para analisar',
        'no_sequence_data': 'Nenhum dado de sequência disponível para análise',
    },
    'en': {
        'title': 'NCBI Query - DNA and Genome',
        'email_label': 'Your Email (required):',
        'email_placeholder': 'your.email@example.com',
        'email_required': 'Please enter your email before searching!',
        'email_invalid': 'Please enter a valid email!',
        'search_label': 'Search Organism:',
        'search_button': '🔍  Search',
        'clear_button': '🗑  Clear',
        'language_label': 'Language:',
        'ready': 'Ready to search',
        'results_cleared': 'Results cleared',
        'searching': '🔍 Searching information about',
        'completed': '✅ Search completed for',
        'error': '❌ Search error:',
        'attention': 'Attention',
        'enter_organism': 'Enter an organism name to search!',
        'error_title': 'Error',
        'error_searching': 'Error searching data:\n',
        'tab_info': '📋  General Information',
        'tab_sequences': '🧬  Sequences',
        'tab_taxonomy': '🌳  Taxonomy',
        'tab_genome': '🔬  Genome',
        'tab_publications': '📚  Publications',
        'searching_nucleotides': 'SEARCHING NUCLEOTIDE SEQUENCES FOR:',
        'total_sequences': '📊 Total sequences found:',
        'showing_first': '📥 Showing first',
        'results': 'results',
        'register': 'RECORD #',
        'annotations': '🔬 Annotations:',
        'sequence': 'SEQUENCE #',
        'size': 'Size:',
        'view_online': '🌐 View online:',
        'sequence_unavailable': '⚠️ Sequence not available (CON/scaffold record without direct sequence)',
        'features': '🧬 Features (first 20):',
        'no_sequences': 'No sequences found.',
        'view_browser': '🌐 VIEW IN NCBI BROWSER',
        'sequences_link': '• Sequences:',
        'genome_complete': '• Complete Genome:',
        'graphics': '• Graphics:',
        'tax_info': 'TAXONOMIC INFORMATION FOR:',
        'organism': 'ORGANISM #',
        'scientific_name': 'Scientific Name:',
        'common_name': 'Common Name:',
        'other_names': 'Other Names:',
        'genetic_code': '🧬 Genetic Code:',
        'no_tax_info': 'No taxonomic information found.',
        'genome_info': 'GENOME INFORMATION FOR:',
        'total_genomes': '📊 Total genomes found:',
        'genome': 'GENOME #',
        'no_genomes': 'No genomes found.',
        'assemblies': 'GENOME ASSEMBLIES',
        'assembly': 'ASSEMBLY #',
        'view_graphics': '🖼️  View Graphics:',
        'visualization_tools': '🌐 NCBI VISUALIZATION TOOLS',
        'visual_resources': '📊 VISUAL RESOURCES:',
        'tip': '💡 Tip: Click the links above to see maps, graphics and interactive images!',
        'gene_info': 'GENE INFORMATION',
        'total_genes': '📊 Total genes found:',
        'gene': 'GENE #',
        'publications': 'RELATED PUBLICATIONS TO:',
        'total_publications': '📊 Total publications found:',
        'publication': 'PUBLICATION #',
        'authors': 'Authors:',
        'summary': 'Summary:',
        'no_publications': 'No publications found.',
        'name': 'Name:',
        'description': 'Description:',
        'division': 'Division:',
        'rank': 'Rank:',
        'lineage': '🌳 Taxonomic Lineage:',
        'organism_label': 'Organism:',
        'status': 'Status:',
        'access': 'Access:',
        'coverage': 'Coverage:',
        'chromosome': 'Chromosome:',
        'location': 'Location:',
        'title_label': 'Title:',
        'year': 'Year:',
        'gene_id': 'Gene ID:',
        'tab_analysis': '📊  Analysis',
        'tab_compare': '🔬  Compare',
        'export_button': '💾  Export',
        'analyze_button': '📊  Analyze',
        'compare_button': '🔍  Compare',
        'nucleotide_composition': 'NUCLEOTIDE COMPOSITION',
        'gc_content': 'GC Content:',
        'at_content': 'AT Content:',
        'total_bases': 'Total Bases:',
        'sequence_1': 'Sequence 1:',
        'sequence_2': 'Sequence 2:',
        'paste_sequence': 'Paste sequence here...',
        'similarity': 'Similarity:',
        'alignment_score': 'Alignment Score:',
        'identical_bases': 'Identical Bases:',
        'export_success': 'Export completed successfully!',
        'success_title': 'Success',
        'error_export': 'Export error',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Choose format:',
        'select_sequence': 'Select a sequence to analyze',
        'no_sequence_data': 'No sequence data available for analysis',
    },
    'es': {
        'title': 'Consulta NCBI - ADN y Genoma',
        'email_label': 'Su Email (requerido):',
        'email_placeholder': 'su.email@ejemplo.com',
        'email_required': '¡Por favor, ingrese su email antes de buscar!',
        'email_invalid': '¡Por favor, ingrese un email válido!',
        'search_label': 'Buscar Organismo:',
        'search_button': '🔍  Buscar',
        'clear_button': '🗑  Limpiar',
        'language_label': 'Idioma:',
        'ready': 'Listo para buscar',
        'results_cleared': 'Resultados eliminados',
        'searching': '🔍 Buscando información sobre',
        'completed': '✅ Búsqueda completada para',
        'error': '❌ Error en la búsqueda:',
        'attention': 'Atención',
        'enter_organism': '¡Ingrese el nombre de un organismo para buscar!',
        'error_title': 'Error',
        'error_searching': 'Error al buscar datos:\n',
        'tab_info': '📋  Información General',
        'tab_sequences': '🧬  Secuencias',
        'tab_taxonomy': '🌳  Taxonomía',
        'tab_genome': '🔬  Genoma',
        'tab_publications': '📚  Publicaciones',
        'searching_nucleotides': 'BUSCANDO SECUENCIAS DE NUCLEÓTIDOS PARA:',
        'total_sequences': '📊 Total de secuencias encontradas:',
        'showing_first': '📥 Mostrando primeros',
        'results': 'resultados',
        'register': 'REGISTRO #',
        'annotations': '🔬 Anotaciones:',
        'sequence': 'SECUENCIA #',
        'size': 'Tamaño:',
        'view_online': '🌐 Ver en línea:',
        'sequence_unavailable': '⚠️ Secuencia no disponible (registro CON/scaffold sin secuencia directa)',
        'features': '🧬 Características (primeras 20):',
        'no_sequences': 'No se encontraron secuencias.',
        'view_browser': '🌐 VER EN NAVEGADOR NCBI',
        'sequences_link': '• Secuencias:',
        'genome_complete': '• Genoma Completo:',
        'graphics': '• Gráficos:',
        'tax_info': 'INFORMACIÓN TAXONÓMICA PARA:',
        'organism': 'ORGANISMO #',
        'scientific_name': 'Nombre Científico:',
        'common_name': 'Nombre Común:',
        'other_names': 'Otros Nombres:',
        'genetic_code': '🧬 Código Genético:',
        'no_tax_info': 'No se encontró información taxonómica.',
        'genome_info': 'INFORMACIÓN DEL GENOMA PARA:',
        'total_genomes': '📊 Total de genomas encontrados:',
        'genome': 'GENOMA #',
        'no_genomes': 'No se encontraron genomas.',
        'assemblies': 'ENSAMBLAJES DE GENOMA',
        'assembly': 'ENSAMBLAJE #',
        'view_graphics': '🖼️  Ver Gráficos:',
        'visualization_tools': '🌐 HERRAMIENTAS DE VISUALIZACIÓN NCBI',
        'visual_resources': '📊 RECURSOS VISUALES:',
        'tip': '💡 Consejo: ¡Haga clic en los enlaces para ver mapas, gráficos e imágenes interactivas!',
        'gene_info': 'INFORMACIÓN DE GENES',
        'total_genes': '📊 Total de genes encontrados:',
        'gene': 'GEN #',
        'publications': 'PUBLICACIONES RELACIONADAS CON:',
        'total_publications': '📊 Total de publicaciones encontradas:',
        'publication': 'PUBLICACIÓN #',
        'authors': 'Autores:',
        'summary': 'Resumen:',
        'no_publications': 'No se encontraron publicaciones.',
        'name': 'Nombre:',
        'description': 'Descripción:',
        'division': 'División:',
        'rank': 'Rango:',
        'lineage': '🌳 Linaje Taxonómico:',
        'organism_label': 'Organismo:',
        'status': 'Estado:',
        'access': 'Acceso:',
        'coverage': 'Cobertura:',
        'chromosome': 'Cromosoma:',
        'location': 'Ubicación:',
        'title_label': 'Título:',
        'year': 'Año:',
        'gene_id': 'ID de Gen:',
        'tab_analysis': '📊  Análisis',
        'tab_compare': '🔬  Comparar',
        'export_button': '💾  Exportar',
        'analyze_button': '📊  Analizar',
        'compare_button': '🔍  Comparar',
        'nucleotide_composition': 'COMPOSICIÓN DE NUCLEÓTIDOS',
        'gc_content': 'Contenido GC:',
        'at_content': 'Contenido AT:',
        'total_bases': 'Total de Bases:',
        'sequence_1': 'Secuencia 1:',
        'sequence_2': 'Secuencia 2:',
        'paste_sequence': 'Pegue la secuencia aquí...',
        'similarity': 'Similitud:',
        'alignment_score': 'Puntuación de Alineación:',
        'identical_bases': 'Bases Idénticas:',
        'export_success': '¡Exportación realizada con éxito!',
        'success_title': 'Éxito',
        'error_export': 'Error al exportar',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Elija el formato:',
        'select_sequence': 'Seleccione una secuencia para analizar',
        'no_sequence_data': 'No hay datos de secuencia disponibles para análisis',
    },
    'fr': {
        'title': 'Consultation NCBI - ADN et Génome',
        'email_label': 'Votre Email (requis):',
        'email_placeholder': 'votre.email@exemple.com',
        'email_required': 'Veuillez entrer votre email avant de rechercher!',
        'email_invalid': 'Veuillez entrer un email valide!',
        'search_label': 'Rechercher un Organisme:',
        'search_button': '🔍  Rechercher',
        'clear_button': '🗑  Effacer',
        'language_label': 'Langue:',
        'ready': 'Prêt à rechercher',
        'results_cleared': 'Résultats effacés',
        'searching': '🔍 Recherche d\'informations sur',
        'completed': '✅ Recherche terminée pour',
        'error': '❌ Erreur de recherche:',
        'attention': 'Attention',
        'enter_organism': 'Entrez le nom d\'un organisme à rechercher!',
        'error_title': 'Erreur',
        'error_searching': 'Erreur lors de la recherche de données:\n',
        'tab_info': '📋  Informations Générales',
        'tab_sequences': '🧬  Séquences',
        'tab_taxonomy': '🌳  Taxonomie',
        'tab_genome': '🔬  Génome',
        'tab_publications': '📚  Publications',
        'searching_nucleotides': 'RECHERCHE DE SÉQUENCES DE NUCLÉOTIDES POUR:',
        'total_sequences': '📊 Total de séquences trouvées:',
        'showing_first': '📥 Affichage des premiers',
        'results': 'résultats',
        'register': 'ENREGISTREMENT #',
        'annotations': '🔬 Annotations:',
        'sequence': 'SÉQUENCE #',
        'size': 'Taille:',
        'view_online': '🌐 Voir en ligne:',
        'sequence_unavailable': '⚠️ Séquence non disponible (enregistrement CON/scaffold sans séquence directe)',
        'features': '🧬 Caractéristiques (20 premières):',
        'no_sequences': 'Aucune séquence trouvée.',
        'view_browser': '🌐 VOIR DANS LE NAVIGATEUR NCBI',
        'sequences_link': '• Séquences:',
        'genome_complete': '• Génome Complet:',
        'graphics': '• Graphiques:',
        'tax_info': 'INFORMATIONS TAXONOMIQUES POUR:',
        'organism': 'ORGANISME #',
        'scientific_name': 'Nom Scientifique:',
        'common_name': 'Nom Commun:',
        'other_names': 'Autres Noms:',
        'genetic_code': '🧬 Code Génétique:',
        'no_tax_info': 'Aucune information taxonomique trouvée.',
        'genome_info': 'INFORMATIONS SUR LE GÉNOME POUR:',
        'total_genomes': '📊 Total de génomes trouvés:',
        'genome': 'GÉNOME #',
        'no_genomes': 'Aucun génome trouvé.',
        'assemblies': 'ASSEMBLAGES DE GÉNOME',
        'assembly': 'ASSEMBLAGE #',
        'view_graphics': '🖼️  Voir les Graphiques:',
        'visualization_tools': '🌐 OUTILS DE VISUALISATION NCBI',
        'visual_resources': '📊 RESSOURCES VISUELLES:',
        'tip': '💡 Astuce: Cliquez sur les liens pour voir des cartes, graphiques et images interactives!',
        'gene_info': 'INFORMATIONS SUR LES GÈNES',
        'total_genes': '📊 Total de gènes trouvés:',
        'gene': 'GÈNE #',
        'publications': 'PUBLICATIONS LIÉES À:',
        'total_publications': '📊 Total de publications trouvées:',
        'publication': 'PUBLICATION #',
        'authors': 'Auteurs:',
        'summary': 'Résumé:',
        'no_publications': 'Aucune publication trouvée.',
        'name': 'Nom:',
        'description': 'Description:',
        'division': 'Division:',
        'rank': 'Rang:',
        'lineage': '🌳 Lignée Taxonomique:',
        'organism_label': 'Organisme:',
        'status': 'Statut:',
        'access': 'Accès:',
        'coverage': 'Couverture:',
        'chromosome': 'Chromosome:',
        'location': 'Emplacement:',
        'title_label': 'Titre:',
        'year': 'Année:',
        'gene_id': 'ID de Gène:',
        'tab_analysis': '📊  Analyse',
        'tab_compare': '🔬  Comparer',
        'export_button': '💾  Exporter',
        'analyze_button': '📊  Analyser',
        'compare_button': '🔍  Comparer',
        'nucleotide_composition': 'COMPOSITION DE NUCLÉOTIDES',
        'gc_content': 'Contenu GC:',
        'at_content': 'Contenu AT:',
        'total_bases': 'Total de Bases:',
        'sequence_1': 'Séquence 1:',
        'sequence_2': 'Séquence 2:',
        'paste_sequence': 'Collez la séquence ici...',
        'similarity': 'Similarité:',
        'alignment_score': 'Score d\'Alignement:',
        'identical_bases': 'Bases Identiques:',
        'export_success': 'Exportation réussie!',
        'success_title': 'Succès',
        'error_export': 'Erreur d\'exportation',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Choisissez le format:',
        'select_sequence': 'Sélectionnez une séquence à analyser',
        'no_sequence_data': 'Aucune donnée de séquence disponible pour l\'analyse',
    },
    'de': {
        'title': 'NCBI-Abfrage - DNA und Genom',
        'email_label': 'Ihre E-Mail (erforderlich):',
        'email_placeholder': 'ihre.email@beispiel.com',
        'email_required': 'Bitte geben Sie Ihre E-Mail ein, bevor Sie suchen!',
        'email_invalid': 'Bitte geben Sie eine gültige E-Mail ein!',
        'search_label': 'Organismus Suchen:',
        'search_button': '🔍  Suchen',
        'clear_button': '🗑  Löschen',
        'language_label': 'Sprache:',
        'ready': 'Bereit zum Suchen',
        'results_cleared': 'Ergebnisse gelöscht',
        'searching': '🔍 Suche nach Informationen über',
        'completed': '✅ Suche abgeschlossen für',
        'error': '❌ Suchfehler:',
        'attention': 'Achtung',
        'enter_organism': 'Geben Sie einen Organismennamen zum Suchen ein!',
        'error_title': 'Fehler',
        'error_searching': 'Fehler bei der Datensuche:\n',
        'tab_info': '📋  Allgemeine Informationen',
        'tab_sequences': '🧬  Sequenzen',
        'tab_taxonomy': '🌳  Taxonomie',
        'tab_genome': '🔬  Genom',
        'tab_publications': '📚  Veröffentlichungen',
        'searching_nucleotides': 'SUCHE NACH NUKLEOTIDSEQUENZEN FÜR:',
        'total_sequences': '📊 Gesamt gefundene Sequenzen:',
        'showing_first': '📥 Zeige erste',
        'results': 'Ergebnisse',
        'register': 'DATENSATZ #',
        'annotations': '🔬 Anmerkungen:',
        'sequence': 'SEQUENZ #',
        'size': 'Größe:',
        'view_online': '🌐 Online ansehen:',
        'sequence_unavailable': '⚠️ Sequenz nicht verfügbar (CON/Scaffold-Eintrag ohne direkte Sequenz)',
        'features': '🧬 Merkmale (erste 20):',
        'no_sequences': 'Keine Sequenzen gefunden.',
        'view_browser': '🌐 IM NCBI-BROWSER ANZEIGEN',
        'sequences_link': '• Sequenzen:',
        'genome_complete': '• Vollständiges Genom:',
        'graphics': '• Grafiken:',
        'tax_info': 'TAXONOMISCHE INFORMATIONEN FÜR:',
        'organism': 'ORGANISMUS #',
        'scientific_name': 'Wissenschaftlicher Name:',
        'common_name': 'Gemeinsamer Name:',
        'other_names': 'Andere Namen:',
        'genetic_code': '🧬 Genetischer Code:',
        'no_tax_info': 'Keine taxonomischen Informationen gefunden.',
        'genome_info': 'GENOMINFORMATIONEN FÜR:',
        'total_genomes': '📊 Gesamt gefundene Genome:',
        'genome': 'GENOM #',
        'no_genomes': 'Keine Genome gefunden.',
        'assemblies': 'GENOMASSEMBLIERUNGEN',
        'assembly': 'ASSEMBLIERUNG #',
        'view_graphics': '🖼️  Grafiken ansehen:',
        'visualization_tools': '🌐 NCBI-VISUALISIERUNGSTOOLS',
        'visual_resources': '📊 VISUELLE RESSOURCEN:',
        'tip': '💡 Tipp: Klicken Sie auf die Links, um Karten, Grafiken und interaktive Bilder zu sehen!',
        'gene_info': 'GENINFORMATIONEN',
        'total_genes': '📊 Gesamt gefundene Gene:',
        'gene': 'GEN #',
        'publications': 'VERWANDTE VERÖFFENTLICHUNGEN ZU:',
        'total_publications': '📊 Gesamt gefundene Veröffentlichungen:',
        'publication': 'VERÖFFENTLICHUNG #',
        'authors': 'Autoren:',
        'summary': 'Zusammenfassung:',
        'no_publications': 'Keine Veröffentlichungen gefunden.',
        'name': 'Name:',
        'description': 'Beschreibung:',
        'division': 'Division:',
        'rank': 'Rang:',
        'lineage': '🌳 Taxonomische Abstammung:',
        'organism_label': 'Organismus:',
        'status': 'Status:',
        'access': 'Zugriff:',
        'coverage': 'Abdeckung:',
        'chromosome': 'Chromosom:',
        'location': 'Standort:',
        'title_label': 'Titel:',
        'year': 'Jahr:',
        'gene_id': 'Gen-ID:',
        'tab_analysis': '📊  Analyse',
        'tab_compare': '🔬  Vergleichen',
        'export_button': '💾  Exportieren',
        'analyze_button': '📊  Analysieren',
        'compare_button': '🔍  Vergleichen',
        'nucleotide_composition': 'NUKLEOTIDZUSAMMENSETZUNG',
        'gc_content': 'GC-Gehalt:',
        'at_content': 'AT-Gehalt:',
        'total_bases': 'Gesamtbasen:',
        'sequence_1': 'Sequenz 1:',
        'sequence_2': 'Sequenz 2:',
        'paste_sequence': 'Sequenz hier einfügen...',
        'similarity': 'Ähnlichkeit:',
        'alignment_score': 'Ausrichtungspunktzahl:',
        'identical_bases': 'Identische Basen:',
        'export_success': 'Export erfolgreich abgeschlossen!',
        'success_title': 'Erfolg',
        'error_export': 'Exportfehler',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Format wählen:',
        'select_sequence': 'Wählen Sie eine Sequenz zum Analysieren',
        'no_sequence_data': 'Keine Sequenzdaten für Analyse verfügbar',
    },
    'zh': {
        'title': 'NCBI查询 - DNA和基因组',
        'email_label': '您的电子邮件（必填）:',
        'email_placeholder': 'your.email@example.com',
        'email_required': '请在搜索前输入您的电子邮件！',
        'email_invalid': '请输入有效的电子邮件！',
        'search_label': '搜索生物体:',
        'search_button': '🔍  搜索',
        'clear_button': '🗑  清除',
        'language_label': '语言:',
        'ready': '准备搜索',
        'results_cleared': '结果已清除',
        'searching': '🔍 正在搜索有关信息',
        'completed': '✅ 搜索完成',
        'error': '❌ 搜索错误:',
        'attention': '注意',
        'enter_organism': '请输入生物体名称进行搜索！',
        'error_title': '错误',
        'error_searching': '搜索数据时出错:\n',
        'tab_info': '📋  一般信息',
        'tab_sequences': '🧬  序列',
        'tab_taxonomy': '🌳  分类学',
        'tab_genome': '🔬  基因组',
        'tab_publications': '📚  出版物',
        'searching_nucleotides': '搜索核苷酸序列:',
        'total_sequences': '📊 找到的序列总数:',
        'showing_first': '📥 显示前',
        'results': '个结果',
        'register': '记录 #',
        'annotations': '🔬 注释:',
        'sequence': '序列 #',
        'size': '大小:',
        'view_online': '🌐 在线查看:',
        'sequence_unavailable': '⚠️ 序列不可用（无直接序列的CON/scaffold记录）',
        'features': '🧬 特征（前20个）:',
        'no_sequences': '未找到序列。',
        'view_browser': '🌐 在NCBI浏览器中查看',
        'sequences_link': '• 序列:',
        'genome_complete': '• 完整基因组:',
        'graphics': '• 图形:',
        'tax_info': '分类信息:',
        'organism': '生物体 #',
        'scientific_name': '学名:',
        'common_name': '俗名:',
        'other_names': '其他名称:',
        'genetic_code': '🧬 遗传密码:',
        'no_tax_info': '未找到分类信息。',
        'genome_info': '基因组信息:',
        'total_genomes': '📊 找到的基因组总数:',
        'genome': '基因组 #',
        'no_genomes': '未找到基因组。',
        'assemblies': '基因组组装',
        'assembly': '组装 #',
        'view_graphics': '🖼️  查看图形:',
        'visualization_tools': '🌐 NCBI可视化工具',
        'visual_resources': '📊 可视化资源:',
        'tip': '💡 提示：点击上面的链接查看地图、图形和交互式图像！',
        'gene_info': '基因信息',
        'total_genes': '📊 找到的基因总数:',
        'gene': '基因 #',
        'publications': '相关出版物:',
        'total_publications': '📊 找到的出版物总数:',
        'publication': '出版物 #',
        'authors': '作者:',
        'summary': '摘要:',
        'no_publications': '未找到出版物。',
        'name': '名称:',
        'description': '描述:',
        'division': '部门:',
        'rank': '等级:',
        'lineage': '🌳 分类谱系:',
        'organism_label': '生物体:',
        'status': '状态:',
        'access': '访问:',
        'coverage': '覆盖率:',
        'chromosome': '染色体:',
        'location': '位置:',
        'title_label': '标题:',
        'year': '年份:',
        'gene_id': '基因ID:',
        'tab_analysis': '📊  分析',
        'tab_compare': '🔬  比较',
        'export_button': '💾  导出',
        'analyze_button': '📊  分析',
        'compare_button': '🔍  比较',
        'nucleotide_composition': '核苷酸组成',
        'gc_content': 'GC含量:',
        'at_content': 'AT含量:',
        'total_bases': '碱基总数:',
        'sequence_1': '序列1:',
        'sequence_2': '序列2:',
        'paste_sequence': '在此粘贴序列...',
        'similarity': '相似度:',
        'alignment_score': '比对得分:',
        'identical_bases': '相同碱基:',
        'export_success': '导出成功完成！',
        'success_title': '成功',
        'error_export': '导出错误',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': '选择格式:',
        'select_sequence': '选择要分析的序列',
        'no_sequence_data': '没有可用于分析的序列数据',
    },
    'ru': {
        'title': 'Запрос NCBI - ДНК и Геном',
        'email_label': 'Ваш Email (обязательно):',
        'email_placeholder': 'your.email@example.com',
        'email_required': 'Пожалуйста, введите свой email перед поиском!',
        'email_invalid': 'Пожалуйста, введите действительный email!',
        'search_label': 'Поиск Организма:',
        'search_button': '🔍  Поиск',
        'clear_button': '🗑  Очистить',
        'language_label': 'Язык:',
        'ready': 'Готов к поиску',
        'results_cleared': 'Результаты очищены',
        'searching': '🔍 Поиск информации о',
        'completed': '✅ Поиск завершен для',
        'error': '❌ Ошибка поиска:',
        'attention': 'Внимание',
        'enter_organism': 'Введите название организма для поиска!',
        'error_title': 'Ошибка',
        'error_searching': 'Ошибка при поиске данных:\n',
        'tab_info': '📋  Общая Информация',
        'tab_sequences': '🧬  Последовательности',
        'tab_taxonomy': '🌳  Таксономия',
        'tab_genome': '🔬  Геном',
        'tab_publications': '📚  Публикации',
        'searching_nucleotides': 'ПОИСК НУКЛЕОТИДНЫХ ПОСЛЕДОВАТЕЛЬНОСТЕЙ ДЛЯ:',
        'total_sequences': '📊 Всего найдено последовательностей:',
        'showing_first': '📥 Показаны первые',
        'results': 'результатов',
        'register': 'ЗАПИСЬ #',
        'annotations': '🔬 Аннотации:',
        'sequence': 'ПОСЛЕДОВАТЕЛЬНОСТЬ #',
        'size': 'Размер:',
        'view_online': '🌐 Просмотреть онлайн:',
        'sequence_unavailable': '⚠️ Последовательность недоступна (запись CON/scaffold без прямой последовательности)',
        'features': '🧬 Особенности (первые 20):',
        'no_sequences': 'Последовательности не найдены.',
        'view_browser': '🌐 ПРОСМОТРЕТЬ В БРАУЗЕРЕ NCBI',
        'sequences_link': '• Последовательности:',
        'genome_complete': '• Полный Геном:',
        'graphics': '• Графика:',
        'tax_info': 'ТАКСОНОМИЧЕСКАЯ ИНФОРМАЦИЯ ДЛЯ:',
        'organism': 'ОРГАНИЗМ #',
        'scientific_name': 'Научное Название:',
        'common_name': 'Общее Название:',
        'other_names': 'Другие Названия:',
        'genetic_code': '🧬 Генетический Код:',
        'no_tax_info': 'Таксономическая информация не найдена.',
        'genome_info': 'ИНФОРМАЦИЯ О ГЕНОМЕ ДЛЯ:',
        'total_genomes': '📊 Всего найдено геномов:',
        'genome': 'ГЕНОМ #',
        'no_genomes': 'Геномы не найдены.',
        'assemblies': 'СБОРКИ ГЕНОМА',
        'assembly': 'СБОРКА #',
        'view_graphics': '🖼️  Просмотреть Графики:',
        'visualization_tools': '🌐 ИНСТРУМЕНТЫ ВИЗУАЛИЗАЦИИ NCBI',
        'visual_resources': '📊 ВИЗУАЛЬНЫЕ РЕСУРСЫ:',
        'tip': '💡 Совет: Нажмите на ссылки выше, чтобы увидеть карты, графики и интерактивные изображения!',
        'gene_info': 'ИНФОРМАЦИЯ О ГЕНАХ',
        'total_genes': '📊 Всего найдено генов:',
        'gene': 'ГЕН #',
        'publications': 'СВЯЗАННЫЕ ПУБЛИКАЦИИ ДЛЯ:',
        'total_publications': '📊 Всего найдено публикаций:',
        'publication': 'ПУБЛИКАЦИЯ #',
        'authors': 'Авторы:',
        'summary': 'Резюме:',
        'no_publications': 'Публикации не найдены.',
        'name': 'Имя:',
        'description': 'Описание:',
        'division': 'Подразделение:',
        'rank': 'Ранг:',
        'lineage': '🌳 Таксономическая Линия:',
        'organism_label': 'Организм:',
        'status': 'Статус:',
        'access': 'Доступ:',
        'coverage': 'Покрытие:',
        'chromosome': 'Хромосома:',
        'location': 'Расположение:',
        'title_label': 'Заголовок:',
        'year': 'Год:',
        'gene_id': 'ID Гена:',
        'tab_analysis': '📊  Анализ',
        'tab_compare': '🔬  Сравнить',
        'export_button': '💾  Экспорт',
        'analyze_button': '📊  Анализировать',
        'compare_button': '🔍  Сравнить',
        'nucleotide_composition': 'СОСТАВ НУКЛЕОТИДОВ',
        'gc_content': 'Содержание GC:',
        'at_content': 'Содержание AT:',
        'total_bases': 'Всего оснований:',
        'sequence_1': 'Последовательность 1:',
        'sequence_2': 'Последовательность 2:',
        'paste_sequence': 'Вставьте последовательность здесь...',
        'similarity': 'Сходство:',
        'alignment_score': 'Оценка выравнивания:',
        'identical_bases': 'Идентичные основания:',
        'export_success': 'Экспорт успешно завершен!',
        'success_title': 'Успех',
        'error_export': 'Ошибка экспорта',
        'export_fasta': 'FASTA',
        'export_pdf': 'PDF',
        'choose_format': 'Выберите формат:',
        'select_sequence': 'Выберите последовательность для анализа',
        'no_sequence_data': 'Нет данных последовательности для анализа',
    },
}

class NCBISearchApp:
    def __init__(self, root):
        self.root = root
        self.current_language = 'pt'  # Idioma padrão
        self.root.title(self.t('title'))
        self.root.geometry("1200x800")
        self.root.configure(bg="#1a1a1a")
        
        # Armazena sequências e dados para análise
        self.sequences = []
        self.current_organism = ""
        self.search_results = {}
        
        # Estilo moderno preto/branco/vermelho
        style = ttk.Style()
        style.theme_use('clam')
        style.configure("TButton", padding=8, relief="flat", background="#c62828", foreground="white", font=("Arial", 10, "bold"))
        style.map("TButton", background=[('active', '#b71c1c')])
        style.configure("TLabel", background="#1a1a1a", foreground="white", font=("Arial", 10))
        style.configure("TNotebook", background="#1a1a1a", borderwidth=0)
        style.configure("TNotebook.Tab", background="#2d2d2d", foreground="white", padding=[20, 10], font=("Arial", 10, "bold"))
        style.map("TNotebook.Tab", 
                 background=[('selected', '#c62828'), ('active', '#c62828')], 
                 foreground=[('selected', 'white'), ('active', 'white')],
                 padding=[('selected', [20, 10]), ('active', [20, 10])])
        
        self.create_widgets()
    
    def t(self, key):
        """Retorna tradução para o idioma atual"""
        return TRANSLATIONS.get(self.current_language, TRANSLATIONS['pt']).get(key, key)
    
    def change_language(self, lang_code):
        """Muda o idioma da interface"""
        self.current_language = lang_code
        self.root.title(self.t('title'))
        # Atualiza labels
        self.email_label.config(text=self.t('email_label'))
        # Atualiza placeholder do email se estiver vazio ou com placeholder
        current_email = self.email_entry.get().strip()
        if not current_email or '@' not in current_email or current_email.endswith('.com') and len(current_email) < 15:
            self.email_entry.delete(0, tk.END)
            self.email_entry.insert(0, self.t('email_placeholder'))
            self.email_entry.config(fg="#999")
        self.search_label.config(text=self.t('search_label'))
        self.search_button.config(text=self.t('search_button'))
        self.clear_button.config(text=self.t('clear_button'))
        self.lang_label.config(text=self.t('language_label'))
        self.status_label.config(text=self.t('ready'))
        # Atualiza abas
        self.notebook.tab(0, text=self.t('tab_info'))
        self.notebook.tab(1, text=self.t('tab_sequences'))
        self.notebook.tab(2, text=self.t('tab_taxonomy'))
        self.notebook.tab(3, text=self.t('tab_genome'))
        self.notebook.tab(4, text=self.t('tab_publications'))
        self.notebook.tab(5, text=self.t('tab_analysis'))
        self.notebook.tab(6, text=self.t('tab_compare'))
        # Atualiza botões
        self.export_button.config(text=self.t('export_button'))
        
    def create_widgets(self):
        # Frame superior - Email e Busca
        top_frame = tk.Frame(self.root, bg="#2d2d2d", padx=15, pady=10)
        top_frame.pack(fill=tk.X, padx=10, pady=10)
        
        # Email
        email_frame = tk.Frame(top_frame, bg="#2d2d2d")
        email_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.email_label = tk.Label(email_frame, text=self.t('email_label'), bg="#2d2d2d", 
                fg="white", font=("Arial", 11, "bold"), width=25, anchor="e")
        self.email_label.pack(side=tk.LEFT, padx=5)
        
        self.email_entry = tk.Entry(email_frame, width=40, font=("Arial", 10), 
                                    bg="#3d3d3d", fg="white", insertbackground="white",
                                    relief=tk.FLAT, bd=2, highlightthickness=1, 
                                    highlightbackground="#555", highlightcolor="#1976d2")
        self.email_entry.pack(side=tk.LEFT, padx=10, ipady=5)
        self.email_entry.insert(0, self.t('email_placeholder'))
        self.email_entry.bind('<FocusIn>', lambda e: self.on_email_focus_in())
        self.email_entry.bind('<FocusOut>', lambda e: self.on_email_focus_out())
        
        # Busca
        search_frame = tk.Frame(top_frame, bg="#2d2d2d")
        search_frame.pack(fill=tk.X)
        
        self.search_label = tk.Label(search_frame, text=self.t('search_label'), bg="#2d2d2d", 
                fg="white", font=("Arial", 11, "bold"), width=25, anchor="e")
        self.search_label.pack(side=tk.LEFT, padx=5)
        
        self.search_entry = tk.Entry(search_frame, width=35, font=("Arial", 11), 
                                     bg="#3d3d3d", fg="white", insertbackground="white",
                                     relief=tk.FLAT, bd=2, highlightthickness=1, 
                                     highlightbackground="#555", highlightcolor="#c62828")
        self.search_entry.pack(side=tk.LEFT, padx=10, ipady=5)
        self.search_entry.bind("<Return>", lambda e: self.search_organism())
        
        self.search_button = tk.Button(search_frame, text=self.t('search_button'), command=self.search_organism,
                 bg="#c62828", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8,
                 activebackground="#b71c1c", activeforeground="white")
        self.search_button.pack(side=tk.LEFT, padx=5)
        
        self.clear_button = tk.Button(search_frame, text=self.t('clear_button'), command=self.clear_results,
                 bg="#424242", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8,
                 activebackground="#616161", activeforeground="white")
        self.clear_button.pack(side=tk.LEFT, padx=5)
        
        self.export_button = tk.Button(search_frame, text=self.t('export_button'), command=self.export_data,
                 bg="#1976d2", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8,
                 activebackground="#1565c0", activeforeground="white")
        self.export_button.pack(side=tk.LEFT, padx=5)
        
        # Seletor de idioma
        self.lang_label = tk.Label(search_frame, text=self.t('language_label'), bg="#2d2d2d",
                                   fg="white", font=("Arial", 10, "bold"))
        self.lang_label.pack(side=tk.LEFT, padx=(20, 5))
        
        languages = [
            ('Português', 'pt'),
            ('English', 'en'),
            ('Español', 'es'),
            ('Français', 'fr'),
            ('Deutsch', 'de'),
            ('中文', 'zh'),
            ('Русский', 'ru')
        ]
        
        self.lang_var = tk.StringVar(value='pt')
        lang_menu = ttk.Combobox(search_frame, textvariable=self.lang_var, values=[lang[0] for lang in languages],
                                state='readonly', width=10, font=("Arial", 9))
        lang_menu.pack(side=tk.LEFT, padx=5)
        lang_menu.bind('<<ComboboxSelected>>', lambda e: self.change_language(languages[lang_menu.current()][1]))
        
        # Frame de banco de dados (removido - simplificar interface)
        
        # Frame de progresso
        self.progress_frame = tk.Frame(self.root, bg="#1a1a1a")
        self.progress_frame.pack(fill=tk.X, padx=15, pady=5)
        
        self.status_label = tk.Label(self.progress_frame, text=self.t('ready'), 
                                     bg="#1a1a1a", fg="#999", font=("Arial", 10, "italic"))
        self.status_label.pack(side=tk.LEFT, padx=5)
        
        self.progress = ttk.Progressbar(self.progress_frame, mode='indeterminate', length=250)
        style = ttk.Style()
        style.configure("TProgressbar", background="#c62828", troughcolor="#2d2d2d", borderwidth=0, thickness=8)
        
        # Notebook (abas) para resultados
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Aba 1: Informações Gerais
        self.info_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.info_frame, text=self.t('tab_info'))
        
        self.info_text = scrolledtext.ScrolledText(self.info_frame, wrap=tk.WORD, 
                                                   font=("Consolas", 10), bg="#0d0d0d", fg="#e0e0e0",
                                                   insertbackground="white", relief=tk.FLAT,
                                                   selectbackground="#c62828", selectforeground="white")
        self.info_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 2: Sequências
        self.seq_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.seq_frame, text=self.t('tab_sequences'))
        
        self.seq_text = scrolledtext.ScrolledText(self.seq_frame, wrap=tk.WORD,
                                                  font=("Consolas", 9), bg="#0d0d0d", fg="#e0e0e0",
                                                  insertbackground="white", relief=tk.FLAT,
                                                  selectbackground="#c62828", selectforeground="white")
        self.seq_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 3: Taxonomia
        self.tax_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.tax_frame, text=self.t('tab_taxonomy'))
        
        self.tax_text = scrolledtext.ScrolledText(self.tax_frame, wrap=tk.WORD,
                                                  font=("Consolas", 10), bg="#0d0d0d", fg="#e0e0e0",
                                                  insertbackground="white", relief=tk.FLAT,
                                                  selectbackground="#c62828", selectforeground="white")
        self.tax_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 4: Genoma
        self.genome_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.genome_frame, text=self.t('tab_genome'))
        
        self.genome_text = scrolledtext.ScrolledText(self.genome_frame, wrap=tk.WORD,
                                                     font=("Consolas", 10), bg="#0d0d0d", fg="#e0e0e0",
                                                     insertbackground="white", relief=tk.FLAT,
                                                     selectbackground="#c62828", selectforeground="white")
        self.genome_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 5: Publicações
        self.pub_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.pub_frame, text=self.t('tab_publications'))
        
        self.pub_text = scrolledtext.ScrolledText(self.pub_frame, wrap=tk.WORD,
                                                  font=("Consolas", 10), bg="#0d0d0d", fg="#e0e0e0",
                                                  insertbackground="white", relief=tk.FLAT,
                                                  selectbackground="#c62828", selectforeground="white")
        self.pub_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 6: Análise
        self.analysis_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.analysis_frame, text=self.t('tab_analysis'))
        
        # Frame para controles de análise
        analysis_control = tk.Frame(self.analysis_frame, bg="#2d2d2d", pady=10)
        analysis_control.pack(fill=tk.X, padx=8, pady=8)
        
        tk.Button(analysis_control, text=self.t('analyze_button'), command=self.analyze_sequences,
                 bg="#c62828", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=20, pady=8).pack(side=tk.LEFT, padx=5)
        
        # Canvas para gráfico
        self.analysis_canvas_frame = tk.Frame(self.analysis_frame, bg="#1a1a1a")
        self.analysis_canvas_frame.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Aba 7: Comparação
        self.compare_frame = tk.Frame(self.notebook, bg="#1a1a1a")
        self.notebook.add(self.compare_frame, text=self.t('tab_compare'))
        
        # Frame superior com inputs
        compare_input_frame = tk.Frame(self.compare_frame, bg="#2d2d2d", pady=10)
        compare_input_frame.pack(fill=tk.X, padx=8, pady=8)
        
        # Sequência 1
        seq1_frame = tk.Frame(compare_input_frame, bg="#2d2d2d")
        seq1_frame.pack(fill=tk.X, pady=5)
        tk.Label(seq1_frame, text=self.t('sequence_1'), bg="#2d2d2d",
                fg="white", font=("Arial", 10, "bold")).pack(side=tk.LEFT, padx=5)
        self.seq1_entry = tk.Text(seq1_frame, height=3, width=80, font=("Consolas", 9),
                                 bg="#3d3d3d", fg="#999", insertbackground="white", relief=tk.FLAT)
        self.seq1_entry.pack(side=tk.LEFT, padx=5)
        self.seq1_entry.insert(1.0, self.t('paste_sequence'))
        self.seq1_entry.bind('<FocusIn>', lambda e: self.on_seq1_focus_in())
        self.seq1_entry.bind('<FocusOut>', lambda e: self.on_seq1_focus_out())
        
        # Sequência 2
        seq2_frame = tk.Frame(compare_input_frame, bg="#2d2d2d")
        seq2_frame.pack(fill=tk.X, pady=5)
        tk.Label(seq2_frame, text=self.t('sequence_2'), bg="#2d2d2d",
                fg="white", font=("Arial", 10, "bold")).pack(side=tk.LEFT, padx=5)
        self.seq2_entry = tk.Text(seq2_frame, height=3, width=80, font=("Consolas", 9),
                                 bg="#3d3d3d", fg="#999", insertbackground="white", relief=tk.FLAT)
        self.seq2_entry.pack(side=tk.LEFT, padx=5)
        self.seq2_entry.insert(1.0, self.t('paste_sequence'))
        self.seq2_entry.bind('<FocusIn>', lambda e: self.on_seq2_focus_in())
        self.seq2_entry.bind('<FocusOut>', lambda e: self.on_seq2_focus_out())
        
        # Frame para botões de alinhamento
        align_buttons_frame = tk.Frame(compare_input_frame, bg="#2d2d2d")
        align_buttons_frame.pack(pady=10)
        
        tk.Label(align_buttons_frame, text="Tipo de alinhamento:", bg="#2d2d2d",
                fg="white", font=("Arial", 10, "bold")).pack(side=tk.LEFT, padx=5)
        
        tk.Button(align_buttons_frame, text="DNA → Proteína (VectorBuilder)", command=lambda: self.compare_sequences('protein'),
                 bg="#4caf50", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8).pack(side=tk.LEFT, padx=5)
        
        tk.Button(align_buttons_frame, text="DNA Local", command=lambda: self.compare_sequences('local'),
                 bg="#1976d2", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8).pack(side=tk.LEFT, padx=5)
        
        tk.Button(align_buttons_frame, text="DNA Global", command=lambda: self.compare_sequences('global'),
                 bg="#c62828", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=15, pady=8).pack(side=tk.LEFT, padx=5)
        
        # Área de resultados
        self.compare_text = scrolledtext.ScrolledText(self.compare_frame, wrap=tk.WORD,
                                                     font=("Consolas", 9), bg="#0d0d0d", fg="#e0e0e0",
                                                     insertbackground="white", relief=tk.FLAT,
                                                     selectbackground="#c62828", selectforeground="white")
        self.compare_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        
        # Torna todos os campos somente leitura
        self.disable_text_editing()
        
    def on_email_focus_in(self):
        """Remove placeholder quando o usuário clica no campo de email"""
        if self.email_entry.get() == self.t('email_placeholder'):
            self.email_entry.delete(0, tk.END)
            self.email_entry.config(fg="white")
    
    def on_email_focus_out(self):
        """Adiciona placeholder se o campo estiver vazio"""
        if not self.email_entry.get().strip():
            self.email_entry.insert(0, self.t('email_placeholder'))
            self.email_entry.config(fg="#999")
    
    def on_seq1_focus_in(self):
        """Remove placeholder quando o usuário clica no campo de sequência 1"""
        current_text = self.seq1_entry.get(1.0, tk.END).strip()
        if current_text == self.t('paste_sequence'):
            self.seq1_entry.delete(1.0, tk.END)
            self.seq1_entry.config(fg="white")
    
    def on_seq1_focus_out(self):
        """Adiciona placeholder se o campo estiver vazio"""
        if not self.seq1_entry.get(1.0, tk.END).strip():
            self.seq1_entry.insert(1.0, self.t('paste_sequence'))
            self.seq1_entry.config(fg="#999")
    
    def on_seq2_focus_in(self):
        """Remove placeholder quando o usuário clica no campo de sequência 2"""
        current_text = self.seq2_entry.get(1.0, tk.END).strip()
        if current_text == self.t('paste_sequence'):
            self.seq2_entry.delete(1.0, tk.END)
            self.seq2_entry.config(fg="white")
    
    def on_seq2_focus_out(self):
        """Adiciona placeholder se o campo estiver vazio"""
        if not self.seq2_entry.get(1.0, tk.END).strip():
            self.seq2_entry.insert(1.0, self.t('paste_sequence'))
            self.seq2_entry.config(fg="#999")
    
    def validate_email(self, email):
        """Valida formato básico de email"""
        import re
        pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        return re.match(pattern, email) is not None
    
    def update_status(self, message, show_progress=False):
        self.status_label.config(text=message)
        if show_progress:
            self.progress.pack(side=tk.LEFT, padx=5)
            self.progress.start(10)
        else:
            self.progress.stop()
            self.progress.pack_forget()
        self.root.update()
        
    def enable_text_editing(self):
        """Habilita edição temporária dos campos de texto"""
        self.info_text.config(state=tk.NORMAL)
        self.seq_text.config(state=tk.NORMAL)
        self.tax_text.config(state=tk.NORMAL)
        self.genome_text.config(state=tk.NORMAL)
        self.pub_text.config(state=tk.NORMAL)
        self.compare_text.config(state=tk.NORMAL)
    
    def disable_text_editing(self):
        """Desabilita edição dos campos de texto"""
        self.info_text.config(state=tk.DISABLED)
        self.seq_text.config(state=tk.DISABLED)
        self.tax_text.config(state=tk.DISABLED)
        self.genome_text.config(state=tk.DISABLED)
        self.pub_text.config(state=tk.DISABLED)
        self.compare_text.config(state=tk.DISABLED)
        
    def clear_results(self):
        """Limpa todos os resultados"""
        self.enable_text_editing()
        self.info_text.delete(1.0, tk.END)
        self.seq_text.delete(1.0, tk.END)
        self.tax_text.delete(1.0, tk.END)
        self.genome_text.delete(1.0, tk.END)
        self.pub_text.delete(1.0, tk.END)
        self.compare_text.delete(1.0, tk.END)
        self.disable_text_editing()
        self.search_entry.delete(0, tk.END)
        self.sequences = []
        self.current_organism = ""
        self.update_status(self.t('results_cleared'))
        
    def search_organism(self):
        """Realiza a busca no NCBI"""
        # Valida email
        email = self.email_entry.get().strip()
        if not email or email == self.t('email_placeholder'):
            messagebox.showwarning(self.t('attention'), self.t('email_required'))
            self.email_entry.focus()
            return
        
        if not self.validate_email(email):
            messagebox.showwarning(self.t('attention'), self.t('email_invalid'))
            self.email_entry.focus()
            return
        
        # Configura email para Entrez
        Entrez.email = email
        
        query = self.search_entry.get().strip()
        if not query:
            messagebox.showwarning(self.t('attention'), self.t('enter_organism'))
            return
        
        self.clear_results()
        self.update_status(f"{self.t('searching')} '{query}'...", True)
        
        # Executa a busca em thread separada para não travar a interface
        thread = threading.Thread(target=self._perform_search, args=(query,), daemon=True)
        thread.start()
    
    def _perform_search(self, query):
        """Realiza a busca em thread separada"""
        try:
            # Limpa sequências anteriores
            self.sequences = []
            self.current_organism = query
            
            # Busca em múltiplos bancos de dados
            self.search_nucleotides(query)
            self.search_taxonomy(query)
            self.search_genome(query)
            self.search_genes(query)
            self.search_publications(query)
            
            self.root.after(0, lambda: self.update_status(f"{self.t('completed')} '{query}' - {datetime.now().strftime('%H:%M:%S')}", False))
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: self.update_status(f"{self.t('error')} {error_msg}", False))
            self.root.after(0, lambda: messagebox.showerror(self.t('error_title'), f"{self.t('error_searching')}{error_msg}"))
    
    def search_nucleotides(self, query):
        """Busca sequências de nucleotídeos"""
        try:
            self.enable_text_editing()
            self.info_text.insert(tk.END, f"{'='*80}\n")
            self.info_text.insert(tk.END, f"{self.t('searching_nucleotides')} {query}\n")
            self.info_text.insert(tk.END, f"{'='*80}\n\n")
            
            # Busca IDs - FILTRADO por organismo específico
            search_query = f'"{query}"[Organism]'
            handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=20, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            count = record["Count"]
            
            self.info_text.insert(tk.END, f"{self.t('total_sequences')} {count}\n")
            self.info_text.insert(tk.END, f"{self.t('showing_first')} {len(id_list)} {self.t('results')}\n\n")
            
            if id_list:
                # Busca detalhes
                handle = Entrez.efetch(db="nucleotide", id=id_list[:10], rettype="gb", retmode="text")
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                
                for i, rec in enumerate(records, 1):
                    # Armazena sequência para análise
                    if str(rec.seq) and len(str(rec.seq)) > 0:
                        self.sequences.append({
                            'id': rec.id,
                            'name': rec.name,
                            'description': rec.description,
                            'sequence': str(rec.seq),
                            'length': len(rec.seq)
                        })
                    
                    self.info_text.insert(tk.END, f"\n{'─'*80}\n")
                    self.info_text.insert(tk.END, f"{self.t('register')}{i}\n")
                    self.info_text.insert(tk.END, f"{'─'*80}\n")
                    self.info_text.insert(tk.END, f"ID: {rec.id}\n")
                    self.info_text.insert(tk.END, f"{self.t('name')} {rec.name}\n")
                    self.info_text.insert(tk.END, f"{self.t('description')} {rec.description}\n")
                    self.info_text.insert(tk.END, f"{self.t('size')} {len(rec.seq)} bp\n")
                    
                    if rec.annotations:
                        self.info_text.insert(tk.END, f"\n{self.t('annotations')}\n")
                        for key, value in list(rec.annotations.items())[:10]:
                            self.info_text.insert(tk.END, f"  • {key}: {value}\n")
                    
                    # Adiciona sequência na aba de sequências
                    self.seq_text.insert(tk.END, f"\n{'='*80}\n")
                    self.seq_text.insert(tk.END, f"{self.t('sequence')}{i}: {rec.id}\n")
                    self.seq_text.insert(tk.END, f"{rec.description}\n")
                    self.seq_text.insert(tk.END, f"{'='*80}\n")
                    self.seq_text.insert(tk.END, f"{self.t('size')} {len(rec.seq)} bp\n")
                    self.seq_text.insert(tk.END, f"{self.t('view_online')} https://www.ncbi.nlm.nih.gov/nuccore/{rec.id}\n\n")
                    
                    # Mostra sequência formatada
                    sequence = str(rec.seq)
                    if sequence and len(sequence) > 0 and sequence != "Seq('')":
                        for j in range(0, min(len(sequence), 5000), 60):
                            self.seq_text.insert(tk.END, f"{sequence[j:j+60]}\n")
                        
                        if len(sequence) > 5000:
                            self.seq_text.insert(tk.END, f"\n... (sequência truncada - total: {len(sequence)} bp)\n")
                    else:
                        self.seq_text.insert(tk.END, f"\n{self.t('sequence_unavailable')}\n")
                    
                    # Features
                    if rec.features:
                        self.seq_text.insert(tk.END, f"\n{self.t('features')}\n")
                        for feat in rec.features[:20]:
                            self.seq_text.insert(tk.END, f"  • {feat.type} - {feat.location}\n")
                            if feat.qualifiers:
                                for key, val in list(feat.qualifiers.items())[:3]:
                                    self.seq_text.insert(tk.END, f"      {key}: {val}\n")
            else:
                self.info_text.insert(tk.END, f"{self.t('no_sequences')}\n")
            
            # Adiciona links para visualização
            self.info_text.insert(tk.END, f"\n\n{'='*80}\n")
            self.info_text.insert(tk.END, f"{self.t('view_browser')}\n")
            self.info_text.insert(tk.END, f"{'='*80}\n")
            query_encoded = query.replace(' ', '+')
            self.info_text.insert(tk.END, f"{self.t('sequences_link')} https://www.ncbi.nlm.nih.gov/nuccore/?term={query_encoded}\n")
            self.info_text.insert(tk.END, f"{self.t('genome_complete')} https://www.ncbi.nlm.nih.gov/genome/?term={query_encoded}\n")
            self.info_text.insert(tk.END, f"{self.t('graphics')} https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/{query_encoded}\n")
                
        except Exception as e:
            self.info_text.insert(tk.END, f"\n❌ Erro ao buscar nucleotídeos: {str(e)}\n")
        finally:
            self.disable_text_editing()
    
    def search_taxonomy(self, query):
        """Busca informações taxonômicas"""
        try:
            self.enable_text_editing()
            self.tax_text.insert(tk.END, f"{'='*80}\n")
            self.tax_text.insert(tk.END, f"{self.t('tax_info')} {query}\n")
            self.tax_text.insert(tk.END, f"{'='*80}\n\n")
            
            # Busca taxonomia
            handle = Entrez.esearch(db="taxonomy", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if id_list:
                # Busca detalhes
                handle = Entrez.efetch(db="taxonomy", id=id_list, retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                for i, tax_rec in enumerate(records, 1):
                    self.tax_text.insert(tk.END, f"\n{'─'*80}\n")
                    self.tax_text.insert(tk.END, f"{self.t('organism')}{i}\n")
                    self.tax_text.insert(tk.END, f"{'─'*80}\n")
                    self.tax_text.insert(tk.END, f"TaxID: {tax_rec.get('TaxId', 'N/A')}\n")
                    self.tax_text.insert(tk.END, f"{self.t('scientific_name')} {tax_rec.get('ScientificName', 'N/A')}\n")
                    
                    if 'OtherNames' in tax_rec:
                        other = tax_rec['OtherNames']
                        if 'GenbankCommonName' in other:
                            self.tax_text.insert(tk.END, f"{self.t('common_name')} {other['GenbankCommonName']}\n")
                        if 'CommonName' in other and isinstance(other['CommonName'], list):
                            self.tax_text.insert(tk.END, f"{self.t('other_names')} {', '.join(other['CommonName'][:5])}\n")
                    
                    self.tax_text.insert(tk.END, f"{self.t('rank')} {tax_rec.get('Rank', 'N/A')}\n")
                    self.tax_text.insert(tk.END, f"{self.t('division')} {tax_rec.get('Division', 'N/A')}\n")
                    
                    if 'LineageEx' in tax_rec:
                        self.tax_text.insert(tk.END, f"\n{self.t('lineage')}\n")
                        for lineage in tax_rec['LineageEx']:
                            self.tax_text.insert(tk.END, f"  • {lineage['Rank']}: {lineage['ScientificName']} (TaxID: {lineage['TaxId']})\n")
                    
                    if 'GeneticCode' in tax_rec:
                        self.tax_text.insert(tk.END, f"\n{self.t('genetic_code')}\n")
                        self.tax_text.insert(tk.END, f"  ID: {tax_rec['GeneticCode'].get('GCId', 'N/A')}\n")
                        self.tax_text.insert(tk.END, f"  Nome: {tax_rec['GeneticCode'].get('GCName', 'N/A')}\n")
                    
                    # Link para visualização
                    tax_id = tax_rec.get('TaxId', '')
                    if tax_id:
                        self.tax_text.insert(tk.END, f"\n🌐 Ver mais: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}\n")
            else:
                self.tax_text.insert(tk.END, f"{self.t('no_tax_info')}\n")
                
        except Exception as e:
            self.tax_text.insert(tk.END, f"\n❌ Erro ao buscar taxonomia: {str(e)}\n")
        finally:
            self.disable_text_editing()
    
    def search_genome(self, query):
        """Busca informações sobre genoma"""
        try:
            self.enable_text_editing()
            self.genome_text.insert(tk.END, f"{'='*80}\n")
            self.genome_text.insert(tk.END, f"{self.t('genome_info')} {query}\n")
            self.genome_text.insert(tk.END, f"{'='*80}\n\n")
            
            # Busca no banco de genomas
            handle = Entrez.esearch(db="genome", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            count = record["Count"]
            
            self.genome_text.insert(tk.END, f"{self.t('total_genomes')} {count}\n\n")
            
            if id_list:
                handle = Entrez.esummary(db="genome", id=",".join(id_list))
                summaries = Entrez.read(handle)
                handle.close()
                
                for i, summary in enumerate(summaries, 1):
                    self.genome_text.insert(tk.END, f"\n{'─'*80}\n")
                    self.genome_text.insert(tk.END, f"{self.t('genome')}{i}\n")
                    self.genome_text.insert(tk.END, f"{'─'*80}\n")
                    
                    for key, value in summary.items():
                        if key not in ['Id']:
                            self.genome_text.insert(tk.END, f"{key}: {value}\n")
            else:
                self.genome_text.insert(tk.END, f"{self.t('no_genomes')}\n")
            
            # Busca também em assembly
            self.genome_text.insert(tk.END, f"\n\n{'='*80}\n")
            self.genome_text.insert(tk.END, f"{self.t('assemblies')}\n")
            self.genome_text.insert(tk.END, f"{'='*80}\n\n")
            
            handle = Entrez.esearch(db="assembly", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            
            if id_list:
                handle = Entrez.esummary(db="assembly", id=",".join(id_list))
                summaries = Entrez.read(handle)
                handle.close()
                
                for i, summary in enumerate(summaries.get('DocumentSummarySet', {}).get('DocumentSummary', []), 1):
                    self.genome_text.insert(tk.END, f"\n{'─'*80}\n")
                    self.genome_text.insert(tk.END, f"{self.t('assembly')}{i}\n")
                    self.genome_text.insert(tk.END, f"{'─'*80}\n")
                    self.genome_text.insert(tk.END, f"{self.t('organism_label')} {summary.get('Organism', 'N/A')}\n")
                    self.genome_text.insert(tk.END, f"{self.t('name')} {summary.get('AssemblyName', 'N/A')}\n")
                    self.genome_text.insert(tk.END, f"{self.t('status')} {summary.get('AssemblyStatus', 'N/A')}\n")
                    accession = summary.get('AssemblyAccession', 'N/A')
                    self.genome_text.insert(tk.END, f"{self.t('access')} {accession}\n")
                    self.genome_text.insert(tk.END, f"{self.t('coverage')} {summary.get('Coverage', 'N/A')}\n")
                    
                    # Link direto para visualização com imagens
                    if accession != 'N/A':
                        self.genome_text.insert(tk.END, f"\n{self.t('view_graphics')} https://www.ncbi.nlm.nih.gov/assembly/{accession}\n")
            
            # Links gerais de visualização
            self.genome_text.insert(tk.END, f"\n\n{'='*80}\n")
            self.genome_text.insert(tk.END, f"{self.t('visualization_tools')}\n")
            self.genome_text.insert(tk.END, f"{'='*80}\n")
            query_encoded = query.replace(' ', '+')
            self.genome_text.insert(tk.END, f"\n{self.t('visual_resources')}\n")
            self.genome_text.insert(tk.END, f"• Genome Browser: https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/{query_encoded}\n")
            self.genome_text.insert(tk.END, f"• Assemblies (com gráficos): https://www.ncbi.nlm.nih.gov/assembly/?term={query_encoded}\n")
            self.genome_text.insert(tk.END, f"• Genome Data Viewer: https://www.ncbi.nlm.nih.gov/genome/gdv/?org={query_encoded}\n")
            self.genome_text.insert(tk.END, f"• Graphics Overview: https://www.ncbi.nlm.nih.gov/genome/?term={query_encoded}\n")
            self.genome_text.insert(tk.END, f"\n{self.t('tip')}\n")
                    
        except Exception as e:
            self.genome_text.insert(tk.END, f"\n❌ Erro ao buscar genoma: {str(e)}\n")
        finally:
            self.disable_text_editing()
    
    def search_genes(self, query):
        """Busca informações sobre genes"""
        try:
            self.enable_text_editing()
            self.info_text.insert(tk.END, f"\n\n{'='*80}\n")
            self.info_text.insert(tk.END, f"{self.t('gene_info')}\n")
            self.info_text.insert(tk.END, f"{'='*80}\n\n")
            
            # Adiciona filtro de organismo para evitar resultados misturados
            gene_query = f"{query}[Organism]"
            handle = Entrez.esearch(db="gene", term=gene_query, retmax=15)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            count = record["Count"]
            
            self.info_text.insert(tk.END, f"{self.t('total_genes')} {count}\n")
            self.info_text.insert(tk.END, f"{self.t('showing_first')} {len(id_list)} {self.t('results')}\n\n")
            
            if id_list:
                handle = Entrez.esummary(db="gene", id=",".join(id_list[:15]))
                summaries = Entrez.read(handle)
                handle.close()
                
                for i, summary in enumerate(summaries.get('DocumentSummarySet', {}).get('DocumentSummary', []), 1):
                    self.info_text.insert(tk.END, f"\n{'─'*60}\n")
                    self.info_text.insert(tk.END, f"{self.t('gene')}{i}\n")
                    self.info_text.insert(tk.END, f"{'─'*60}\n")
                    self.info_text.insert(tk.END, f"{self.t('gene_id')} {summary.get('Id', 'N/A')}\n")
                    self.info_text.insert(tk.END, f"{self.t('name')} {summary.get('Name', 'N/A')}\n")
                    self.info_text.insert(tk.END, f"{self.t('description')} {summary.get('Description', 'N/A')}\n")
                    self.info_text.insert(tk.END, f"{self.t('organism_label')} {summary.get('Organism', {}).get('ScientificName', 'N/A')}\n")
                    self.info_text.insert(tk.END, f"{self.t('chromosome')} {summary.get('Chromosome', 'N/A')}\n")
                    self.info_text.insert(tk.END, f"{self.t('location')} {summary.get('MapLocation', 'N/A')}\n")
                    
        except Exception as e:
            self.info_text.insert(tk.END, f"\n❌ Erro ao buscar genes: {str(e)}\n")
        finally:
            self.disable_text_editing()
    
    def search_publications(self, query):
        """Busca publicações relacionadas"""
        try:
            self.enable_text_editing()
            self.pub_text.insert(tk.END, f"{'='*80}\n")
            self.pub_text.insert(tk.END, f"{self.t('publications')} {query}\n")
            self.pub_text.insert(tk.END, f"{'='*80}\n\n")
            
            handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            count = record["Count"]
            
            self.pub_text.insert(tk.END, f"{self.t('total_publications')} {count}\n")
            self.pub_text.insert(tk.END, f"{self.t('showing_first')} {len(id_list)} {self.t('results')}\n\n")
            
            if id_list:
                handle = Entrez.efetch(db="pubmed", id=id_list[:20], retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                for i, article in enumerate(records.get('PubmedArticle', []), 1):
                    medline = article.get('MedlineCitation', {})
                    art = medline.get('Article', {})
                    
                    self.pub_text.insert(tk.END, f"\n{'─'*80}\n")
                    self.pub_text.insert(tk.END, f"{self.t('publication')}{i}\n")
                    self.pub_text.insert(tk.END, f"{'─'*80}\n")
                    
                    pmid = medline.get('PMID', 'N/A')
                    self.pub_text.insert(tk.END, f"PMID: {pmid}\n")
                    
                    title = art.get('ArticleTitle', 'N/A')
                    self.pub_text.insert(tk.END, f"{self.t('title_label')} {title}\n")
                    
                    # Autores
                    authors = art.get('AuthorList', [])
                    if authors:
                        author_names = []
                        for author in authors[:5]:
                            last = author.get('LastName', '')
                            init = author.get('Initials', '')
                            if last:
                                author_names.append(f"{last} {init}")
                        if author_names:
                            self.pub_text.insert(tk.END, f"{self.t('authors')} {', '.join(author_names)}")
                            if len(authors) > 5:
                                self.pub_text.insert(tk.END, f" et al. ({len(authors)} autores)")
                            self.pub_text.insert(tk.END, "\n")
                    
                    # Journal
                    journal = art.get('Journal', {})
                    journal_title = journal.get('Title', 'N/A')
                    self.pub_text.insert(tk.END, f"Journal: {journal_title}\n")
                    
                    # Data
                    pub_date = journal.get('JournalIssue', {}).get('PubDate', {})
                    year = pub_date.get('Year', 'N/A')
                    self.pub_text.insert(tk.END, f"{self.t('year')} {year}\n")
                    
                    # Abstract
                    abstract = art.get('Abstract', {}).get('AbstractText', [])
                    if abstract:
                        self.pub_text.insert(tk.END, f"\n{self.t('summary')}\n")
                        if isinstance(abstract, list):
                            for abs_text in abstract[:1]:
                                self.pub_text.insert(tk.END, f"{str(abs_text)[:500]}...\n")
                        else:
                            self.pub_text.insert(tk.END, f"{str(abstract)[:500]}...\n")
                    
                    self.pub_text.insert(tk.END, f"\nLink: https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n")
            else:
                self.pub_text.insert(tk.END, f"{self.t('no_publications')}\n")
                
        except Exception as e:
            self.pub_text.insert(tk.END, f"\n❌ Erro ao buscar publicações: {str(e)}\n")
        finally:
            self.disable_text_editing()


    def analyze_sequences(self):
        """Analisa composição de nucleotídeos e gera gráficos"""
        if not self.sequences:
            messagebox.showinfo(self.t('attention'), self.t('no_sequence_data'))
            return
        
        # Limpa canvas anterior
        for widget in self.analysis_canvas_frame.winfo_children():
            widget.destroy()
        
        # Pega primeira sequência válida
        seq_data = self.sequences[0]
        sequence = seq_data['sequence'].upper()
        
        # Conta nucleotídeos
        counts = Counter(sequence)
        total = len(sequence)
        
        a_count = counts.get('A', 0)
        t_count = counts.get('T', 0)
        g_count = counts.get('G', 0)
        c_count = counts.get('C', 0)
        
        gc_content = ((g_count + c_count) / total * 100) if total > 0 else 0
        at_content = ((a_count + t_count) / total * 100) if total > 0 else 0
        
        # Cria figura com 2 subplots
        fig = Figure(figsize=(12, 5), facecolor='#1a1a1a')
        
        # Gráfico de pizza
        ax1 = fig.add_subplot(121)
        colors_pie = ['#ff6b6b', '#4ecdc4', '#ffe66d', '#a8e6cf']
        sizes = [a_count, t_count, g_count, c_count]
        labels = [f'A: {a_count} ({a_count/total*100:.1f}%)',
                 f'T: {t_count} ({t_count/total*100:.1f}%)',
                 f'G: {g_count} ({g_count/total*100:.1f}%)',
                 f'C: {c_count} ({c_count/total*100:.1f}%)']
        
        ax1.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
               startangle=90, textprops={'color': 'white', 'fontsize': 10})
        ax1.set_title(self.t('nucleotide_composition'), color='white', fontsize=12, pad=20)
        ax1.set_facecolor('#1a1a1a')
        
        # Gráfico de barras
        ax2 = fig.add_subplot(122)
        bases = ['A', 'T', 'G', 'C']
        counts_list = [a_count, t_count, g_count, c_count]
        bars = ax2.bar(bases, counts_list, color=colors_pie, edgecolor='white', linewidth=1.5)
        
        ax2.set_title(f'{self.t("total_bases")} {total}', color='white', fontsize=12, pad=20)
        ax2.set_xlabel('Base', color='white', fontsize=10)
        ax2.set_ylabel('Contagem', color='white', fontsize=10)
        ax2.tick_params(colors='white')
        ax2.set_facecolor('#1a1a1a')
        ax2.spines['bottom'].set_color('white')
        ax2.spines['top'].set_color('white')
        ax2.spines['right'].set_color('white')
        ax2.spines['left'].set_color('white')
        
        # Adiciona valores nas barras
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom', color='white', fontsize=10)
        
        fig.tight_layout()
        
        # Adiciona canvas ao frame
        canvas = FigureCanvasTkAgg(fig, master=self.analysis_canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Adiciona informações textuais
        info_frame = tk.Frame(self.analysis_canvas_frame, bg="#2d2d2d", pady=15)
        info_frame.pack(fill=tk.X, padx=20)
        
        info_text = f"""
        {self.t('gc_content')} {gc_content:.2f}%
        {self.t('at_content')} {at_content:.2f}%
        {self.t('total_bases')} {total}
        
        Sequência: {seq_data['id']} - {seq_data['description'][:60]}...
        """
        
        tk.Label(info_frame, text=info_text, bg="#2d2d2d", fg="white",
                font=("Arial", 11), justify=tk.LEFT).pack(padx=10)
    
    def _process_protein_alignment(self, aligned_seq1, aligned_seq2, score, protein1, protein2,
                                   align_type_label, gap_open, gap_extend):
        """Processa e exibe alinhamento de proteínas"""
        # Calcula estatísticas para proteínas
        matches = 0
        similar = 0  # Substituições conservativas
        gaps_count = 0
        
        # Matriz de similaridade BLOSUM62 simplificada para determinar similaridade
        blosum_similar = {
            'positive': ['RKH', 'DE', 'AVLIM', 'FYW', 'STNQ'],  # Grupos de aminoácidos similares
        }
        
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b and a != '-':
                matches += 1
            elif a == '-' or b == '-':
                gaps_count += 1
            else:
                # Verifica se é substituição conservativa
                is_similar = False
                for group in blosum_similar['positive']:
                    if a in group and b in group:
                        is_similar = True
                        break
                if is_similar:
                    similar += 1
        
        gaps1 = aligned_seq1.count('-')
        gaps2 = aligned_seq2.count('-')
        total_gaps = gaps1 + gaps2
        alignment_length = len(aligned_seq1)
        
        # Identity: apenas matches perfeitos
        identity = (matches / alignment_length) * 100
        
        # Similarity: matches + substituições conservativas
        similarity = ((matches + similar) / alignment_length) * 100
        
        # Salva dados para exportação PDF
        self.last_alignment_data = {
            'type': f'{align_type_label} (BLOSUM62)',
            'score': f'{score}',
            'identity': f'{matches}/{alignment_length} ({identity:.2f}%)',
            'similarity': f'{matches+similar}/{alignment_length} ({similarity:.2f}%)',
            'gaps': f'{total_gaps}/{alignment_length} ({(total_gaps/alignment_length*100):.2f}%)'
        }
        
        # Mostra estatísticas
        self.compare_text.insert(tk.END, "📊 ESTATÍSTICAS DO ALINHAMENTO DE PROTEÍNA:\n")
        self.compare_text.insert(tk.END, f"{'─'*90}\n")
        self.compare_text.insert(tk.END, "  PARÂMETROS:\n")
        self.compare_text.insert(tk.END, f"    Tipo: {align_type_label}\n")
        self.compare_text.insert(tk.END, f"    Matriz: BLOSUM62\n")
        self.compare_text.insert(tk.END, f"    Gap open penalty: {gap_open}\n")
        self.compare_text.insert(tk.END, f"    Gap extend penalty: {gap_extend}\n\n")
        self.compare_text.insert(tk.END, "  RESULTADOS:\n")
        self.compare_text.insert(tk.END, f"    Score de Alinhamento: {score}\n")
        self.compare_text.insert(tk.END, f"    Sequence 1 protein length: {len(protein1)} aa\n")
        self.compare_text.insert(tk.END, f"    Sequence 2 protein length: {len(protein2)} aa\n")
        self.compare_text.insert(tk.END, f"    Alignment length: {alignment_length} aa\n\n")
        self.compare_text.insert(tk.END, "  SIMILARIDADE E IDENTIDADE:\n")
        self.compare_text.insert(tk.END, f"    Similarity:     {matches+similar}/{alignment_length} ({similarity:.2f}%)\n")
        self.compare_text.insert(tk.END, f"    Identity:       {matches}/{alignment_length} ({identity:.2f}%)\n")
        self.compare_text.insert(tk.END, f"    Identical:      {matches}\n")
        self.compare_text.insert(tk.END, f"    Similar:        {similar}\n")
        self.compare_text.insert(tk.END, f"    Gaps:           {total_gaps}/{alignment_length} ({(total_gaps/alignment_length*100):.2f}%) (Seq1: {gaps1}, Seq2: {gaps2})\n")
        self.compare_text.insert(tk.END, f"{'─'*90}\n\n")
        
        self.compare_text.insert(tk.END, "🧬 ALINHAMENTO VISUAL DE PROTEÍNA:\n")
        self.compare_text.insert(tk.END, f"{'─'*90}\n\n")
        
        # Mostra alinhamento em blocos de 60 aminoácidos
        block_size = 60
        for start in range(0, alignment_length, block_size):
            end = min(start + block_size, alignment_length)
            block_seq1 = aligned_seq1[start:end]
            block_seq2 = aligned_seq2[start:end]
            
            # Cria linha de símbolos (simplificado para proteína)
            match_line = ''
            for a, b in zip(block_seq1, block_seq2):
                if a == b and a != '-':
                    match_line += '|'
                elif a == '-' or b == '-':
                    match_line += ' '
                else:
                    # Verifica similaridade
                    is_sim = False
                    for group in blosum_similar['positive']:
                        if a in group and b in group:
                            is_sim = True
                            break
                    match_line += ':' if is_sim else '.'
            
            # Formata em grupos de 10
            def format_with_spaces(seq):
                groups = [seq[i:i+10] for i in range(0, len(seq), 10)]
                return ' '.join(groups)
            
            formatted_seq1 = format_with_spaces(block_seq1)
            formatted_seq2 = format_with_spaces(block_seq2)
            formatted_match = format_with_spaces(match_line)
            
            # Calcula posição real (sem gaps)
            pos1 = sum(1 for c in aligned_seq1[:start] if c != '-') + 1
            pos2 = sum(1 for c in aligned_seq2[:start] if c != '-') + 1
            
            # Mostra o bloco
            self.compare_text.insert(tk.END, f"Prot1: {pos1:>4}  {formatted_seq1}\n")
            self.compare_text.insert(tk.END, f"            {formatted_match}\n")
            self.compare_text.insert(tk.END, f"Prot2: {pos2:>4}  {formatted_seq2}\n\n")
        
        # Legenda
        self.compare_text.insert(tk.END, f"\n{'─'*90}\n")
        self.compare_text.insert(tk.END, "📖 LEGENDA (Alinhamento de Proteína):\n")
        self.compare_text.insert(tk.END, "  |  = Aminoácido idêntico\n")
        self.compare_text.insert(tk.END, "  :  = Aminoácidos similares (substituição conservativa)\n")
        self.compare_text.insert(tk.END, "  .  = Aminoácidos diferentes\n")
        self.compare_text.insert(tk.END, "     = Gap\n\n")
        self.compare_text.insert(tk.END, "  DIFERENÇA:\n")
        self.compare_text.insert(tk.END, "  • Similarity = Idênticos + Similares (como VectorBuilder)\n")
        self.compare_text.insert(tk.END, "  • Identity = Apenas idênticos\n")
        self.compare_text.insert(tk.END, f"{'─'*90}\n")
    
    def compare_sequences(self, alignment_type='local'):
        """Compara duas sequências e mostra alinhamento formatado estilo VectorBuilder"""
        seq1 = self.seq1_entry.get(1.0, tk.END).strip().upper()
        seq2 = self.seq2_entry.get(1.0, tk.END).strip().upper()
        
        # Remove texto placeholder
        if seq1 == self.t('paste_sequence').upper():
            seq1 = ""
        if seq2 == self.t('paste_sequence').upper():
            seq2 = ""
        
        if not seq1 or not seq2:
            messagebox.showwarning(self.t('attention'), 'Por favor, insira ambas as sequências')
            return
        
        # Remove caracteres não-nucleotídeos e espaços
        seq1 = ''.join(c for c in seq1 if c in 'ATGC')
        seq2 = ''.join(c for c in seq2 if c in 'ATGC')
        
        if len(seq1) == 0 or len(seq2) == 0:
            messagebox.showwarning(self.t('attention'), 'Sequências inválidas')
            return
        
        self.compare_text.config(state=tk.NORMAL)
        self.compare_text.delete(1.0, tk.END)
        
        self.compare_text.insert(tk.END, f"{'='*90}\n")
        align_type_label = "PROTEÍNA (DNA traduzido)" if alignment_type == 'protein' else ("LOCAL" if alignment_type == 'local' else "GLOBAL")
        self.compare_text.insert(tk.END, f"COMPARAÇÃO DE SEQUÊNCIAS - ALINHAMENTO {align_type_label}\n")
        self.compare_text.insert(tk.END, f"{'='*90}\n\n")
        
        # Parâmetros de alinhamento (similares ao VectorBuilder)
        match_score = 2      # Score para match
        mismatch_score = -1  # Penalidade para mismatch
        gap_open = -2.0      # Penalidade para abrir gap (VectorBuilder usa -2.0)
        gap_extend = -2.0    # Penalidade para estender gap (VectorBuilder usa -2.0)
        
        # Escolhe o tipo de alinhamento
        if alignment_type == 'protein':
            # Alinhamento baseado em proteína traduzida (como VectorBuilder)
            from Bio.Seq import Seq
            from Bio import Align
            
            self.compare_text.insert(tk.END, "ℹ️  ALINHAMENTO BASEADO EM PROTEÍNA TRADUZIDA (VectorBuilder style):\n")
            self.compare_text.insert(tk.END, "   1. Testa todos os 6 frames de tradução (3 forward + 3 reverse)\n")
            self.compare_text.insert(tk.END, "   2. Escolhe o frame com melhor score de alinhamento\n")
            self.compare_text.insert(tk.END, "   3. Alinha usando matriz BLOSUM62\n\n")
            
            # Testa todos os 6 frames de tradução e escolhe o melhor
            try:
                # Configurar alinhador uma vez
                aligner = Align.PairwiseAligner()
                aligner.mode = 'global'  # Alinhamento global (VectorBuilder style)
                aligner.open_gap_score = gap_open
                aligner.extend_gap_score = gap_extend
                # Gaps nas extremidades não são penalizados (VectorBuilder style)
                aligner.target_end_gap_score = 0.0
                aligner.query_end_gap_score = 0.0
                aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
                
                best_score = float('-inf')
                best_alignment_data = None
                best_frame_info = None
                
                # Testar 6 combinações de frames (3 para seq1 × 2 orientações para seq2)
                for frame1 in range(3):
                    seq1_frame = seq1[frame1:]
                    # Ajusta para múltiplo de 3
                    seq1_frame = seq1_frame[:len(seq1_frame) - (len(seq1_frame) % 3)]
                    protein1 = str(Seq(seq1_frame).translate(to_stop=False))
                    
                    for frame2 in range(3):
                        seq2_frame = seq2[frame2:]
                        seq2_frame = seq2_frame[:len(seq2_frame) - (len(seq2_frame) % 3)]
                        protein2 = str(Seq(seq2_frame).translate(to_stop=False))
                        
                        # Alinha proteínas
                        try:
                            alignment = next(aligner.align(protein1, protein2))
                            if alignment.score > best_score:
                                best_score = alignment.score
                                best_alignment_data = (str(alignment[0]), str(alignment[1]), alignment.score)
                                best_frame_info = (frame1 + 1, frame2 + 1, protein1, protein2, len(seq1_frame), len(seq2_frame))
                        except StopIteration:
                            continue
                    
                    # Também testa reverse complement da seq2
                    seq2_rc = str(Seq(seq2).reverse_complement())
                    for frame2 in range(3):
                        seq2_frame = seq2_rc[frame2:]
                        seq2_frame = seq2_frame[:len(seq2_frame) - (len(seq2_frame) % 3)]
                        protein2 = str(Seq(seq2_frame).translate(to_stop=False))
                        
                        try:
                            alignment = next(aligner.align(protein1, protein2))
                            if alignment.score > best_score:
                                best_score = alignment.score
                                best_alignment_data = (str(alignment[0]), str(alignment[1]), alignment.score)
                                best_frame_info = (frame1 + 1, -(frame2 + 1), protein1, protein2, len(seq1_frame), len(seq2_frame))
                        except StopIteration:
                            continue
                
                if best_alignment_data is None:
                    self.compare_text.insert(tk.END, "❌ Não foi possível criar alinhamento em nenhum frame.\n")
                else:
                    aligned_seq1, aligned_seq2, score = best_alignment_data
                    frame1, frame2, protein1, protein2, bp1, bp2 = best_frame_info
                    
                    strand2 = "forward" if frame2 > 0 else "reverse"
                    self.compare_text.insert(tk.END, f"✓ Melhor alinhamento: Frame {frame1} (Seq1) × Frame {abs(frame2)} (Seq2-{strand2})\n")
                    self.compare_text.insert(tk.END, f"DNA Seq1 ({bp1} bp) → Proteína ({len(protein1)} aa)\n")
                    self.compare_text.insert(tk.END, f"DNA Seq2 ({bp2} bp) → Proteína ({len(protein2)} aa)\n\n")
                    
                    # Processa alinhamento de proteínas (usa função separada)
                    self._process_protein_alignment(aligned_seq1, aligned_seq2, score, protein1, protein2, 
                                                    align_type_label, gap_open, gap_extend)
                
            except Exception as e:
                self.compare_text.insert(tk.END, f"❌ Erro na tradução/alinhamento: {str(e)}\n")
            
        elif alignment_type == 'local':
            # Alinhamento local - encontra regiões de melhor similaridade
            alignments = pairwise2.align.localms(seq1, seq2, 
                                                  match_score, mismatch_score,
                                                  gap_open, gap_extend)
            self.compare_text.insert(tk.END, "ℹ️  ALINHAMENTO LOCAL: Encontra as regiões de maior similaridade.\n")
            self.compare_text.insert(tk.END, "   Ideal para sequências de espécies diferentes ou tamanhos muito distintos.\n\n")
        else:
            # Alinhamento global - alinha de ponta a ponta
            alignments = pairwise2.align.globalms(seq1, seq2, 
                                                   match_score, mismatch_score,
                                                   gap_open, gap_extend)
            self.compare_text.insert(tk.END, "ℹ️  ALINHAMENTO GLOBAL: Alinha as sequências de ponta a ponta.\n")
            self.compare_text.insert(tk.END, "   Ideal para sequências da mesma espécie e tamanhos similares.\n\n")
        
        if alignment_type != 'protein':  # DNA direto
            alignments = pairwise2.align.globalms(seq1, seq2, 
                                                   match_score, mismatch_score,
                                                   gap_open, gap_extend) if alignment_type == 'global' else \
                         pairwise2.align.localms(seq1, seq2, 
                                                  match_score, mismatch_score,
                                                  gap_open, gap_extend)
        
            if alignments:
                best_alignment = alignments[0]
                aligned_seq1, aligned_seq2, score, begin, end = best_alignment
                
                # Calcula estatísticas detalhadas
                matches = 0
                transitions = 0
                transversions = 0
                gaps_count = 0
                
                purines = ['A', 'G']
                pyrimidines = ['C', 'T']
                
                for a, b in zip(aligned_seq1, aligned_seq2):
                    if a == b and a != '-':
                        matches += 1
                    elif a == '-' or b == '-':
                        gaps_count += 1
                    elif (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                        transitions += 1  # Substituição conservativa
                    else:
                        transversions += 1  # Substituição não-conservativa
                
                gaps1 = aligned_seq1.count('-')
                gaps2 = aligned_seq2.count('-')
                total_gaps = gaps1 + gaps2
                alignment_length = len(aligned_seq1)
                
                # Identity: apenas matches perfeitos
                identity = (matches / alignment_length) * 100
                
                # Similarity: matches + transições (substituições conservativas)
                similarity = ((matches + transitions) / alignment_length) * 100
                similarity_original = (matches / max(len(seq1), len(seq2))) * 100
                
                # Mostra estatísticas
                self.compare_text.insert(tk.END, "📊 ESTATÍSTICAS DO ALINHAMENTO:\n")
                self.compare_text.insert(tk.END, f"{'─'*90}\n")
                self.compare_text.insert(tk.END, "  PARÂMETROS:\n")
                self.compare_text.insert(tk.END, f"    Tipo: {align_type_label}\n")
                self.compare_text.insert(tk.END, f"    Match score: {match_score}\n")
                self.compare_text.insert(tk.END, f"    Mismatch score: {mismatch_score}\n")
                self.compare_text.insert(tk.END, f"    Gap open penalty: {gap_open}\n")
                self.compare_text.insert(tk.END, f"    Gap extend penalty: {gap_extend}\n\n")
                self.compare_text.insert(tk.END, "  RESULTADOS:\n")
                self.compare_text.insert(tk.END, f"    {self.t('alignment_score')} {score}\n")
                self.compare_text.insert(tk.END, f"    Sequence 1 length: {len(seq1)} bp\n")
                self.compare_text.insert(tk.END, f"    Sequence 2 length: {len(seq2)} bp\n")
                self.compare_text.insert(tk.END, f"    Alignment length: {alignment_length} bp\n\n")
                self.compare_text.insert(tk.END, "  SIMILARIDADE E IDENTIDADE:\n")
                self.compare_text.insert(tk.END, f"    Similarity:     {matches+transitions}/{alignment_length} ({similarity:.2f}%)\n")
                self.compare_text.insert(tk.END, f"    Identity:       {matches}/{alignment_length} ({identity:.2f}%)\n")
                self.compare_text.insert(tk.END, f"    Matches (|):    {matches}\n")
                self.compare_text.insert(tk.END, f"    Transitions (:): {transitions}\n")
                self.compare_text.insert(tk.END, f"    Transversions (.): {transversions}\n")
                self.compare_text.insert(tk.END, f"    Gaps:           {total_gaps}/{alignment_length} ({(total_gaps/alignment_length*100):.2f}%) (Seq1: {gaps1}, Seq2: {gaps2})\n")
                
                # Salva dados para exportação PDF
                self.last_alignment_data = {
                    'type': f'{align_type_label} (DNA)',
                    'score': f'{score}',
                    'identity': f'{matches}/{alignment_length} ({identity:.2f}%)',
                    'similarity': f'{matches+transitions}/{alignment_length} ({similarity:.2f}%)',
                    'gaps': f'{total_gaps}/{alignment_length} ({(total_gaps/alignment_length*100):.2f}%)'
                }
                
                # Adiciona nota sobre qualidade do alinhamento
                if alignment_type == 'local':
                    coverage1 = (alignment_length - gaps1) / len(seq1) * 100
                    coverage2 = (alignment_length - gaps2) / len(seq2) * 100
                    self.compare_text.insert(tk.END, f"    Cobertura Seq1: {coverage1:.2f}% | Cobertura Seq2: {coverage2:.2f}%\n")
                
                self.compare_text.insert(tk.END, f"{'─'*90}\n\n")
                
                self.compare_text.insert(tk.END, "🧬 ALINHAMENTO VISUAL:\n")
                self.compare_text.insert(tk.END, f"{'─'*90}\n\n")
                
                # Função auxiliar para criar símbolos de correspondência
                def get_match_symbol(a, b):
                    if a == b and a != '-':
                        return '|'  # Match perfeito
                    elif a == '-' or b == '-':
                        return ' '  # Gap
                    else:
                        # Verifica se é transição ou transversão
                        purines = ['A', 'G']
                        pyrimidines = ['C', 'T']
                        if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                            return ':'  # Transição (purina-purina ou pirimidina-pirimidina)
                        else:
                            return '.'  # Transversão
                
                # Mostra alinhamento em blocos de 60 caracteres (6 grupos de 10)
                block_size = 60
                for start in range(0, alignment_length, block_size):
                    end = min(start + block_size, alignment_length)
                    block_seq1 = aligned_seq1[start:end]
                    block_seq2 = aligned_seq2[start:end]
                    
                    # Cria linha de símbolos de correspondência
                    match_line = ''.join(get_match_symbol(a, b) for a, b in zip(block_seq1, block_seq2))
                    
                    # Formata sequências em grupos de 10
                    def format_with_spaces(seq):
                        groups = [seq[i:i+10] for i in range(0, len(seq), 10)]
                        return ' '.join(groups)
                    
                    formatted_seq1 = format_with_spaces(block_seq1)
                    formatted_seq2 = format_with_spaces(block_seq2)
                    formatted_match = format_with_spaces(match_line)
                    
                    # Calcula posição real (sem contar gaps)
                    pos1 = sum(1 for c in aligned_seq1[:start] if c != '-') + 1
                    pos2 = sum(1 for c in aligned_seq2[:start] if c != '-') + 1
                    
                    # Mostra o bloco formatado
                    self.compare_text.insert(tk.END, f"Seq1: {pos1:>5}  {formatted_seq1}\n")
                    self.compare_text.insert(tk.END, f"             {formatted_match}\n")
                    self.compare_text.insert(tk.END, f"Seq2: {pos2:>5}  {formatted_seq2}\n\n")
                
                # Legenda dos símbolos
                self.compare_text.insert(tk.END, f"\n{'─'*90}\n")
                self.compare_text.insert(tk.END, "📖 LEGENDA:\n")
                self.compare_text.insert(tk.END, "  SÍMBOLOS DE ALINHAMENTO:\n")
                self.compare_text.insert(tk.END, "    |  = Match perfeito (bases idênticas)\n")
                self.compare_text.insert(tk.END, "    :  = Transição (A↔G ou C↔T - substituição conservativa)\n")
                self.compare_text.insert(tk.END, "    .  = Transversão (purina↔pirimidina - substituição radical)\n")
                self.compare_text.insert(tk.END, "       = Mismatch ou gap\n\n")
                self.compare_text.insert(tk.END, "  DIFERENÇA ENTRE SIMILARITY E IDENTITY:\n")
                self.compare_text.insert(tk.END, "    • Similarity: Inclui matches (|) + transições (:)\n")
                self.compare_text.insert(tk.END, "      → Considera substituições conservativas como 'similares'\n")
                self.compare_text.insert(tk.END, "    • Identity: Apenas matches perfeitos (|)\n")
                self.compare_text.insert(tk.END, "      → Requer bases exatamente iguais\n")
                self.compare_text.insert(tk.END, f"{'─'*90}\n")
            
        # Fecha o texto como somente leitura
        self.compare_text.config(state=tk.DISABLED)
    
    def export_data(self):
        """Exporta dados em FASTA ou PDF"""
        if not self.sequences:
            messagebox.showwarning(self.t('attention'), self.t('no_sequence_data'))
            return
        
        # Cria janela de diálogo para escolher formato
        export_window = tk.Toplevel(self.root)
        export_window.title(self.t('export_button'))
        export_window.geometry("300x150")
        export_window.configure(bg="#1a1a1a")
        export_window.transient(self.root)
        export_window.grab_set()
        
        tk.Label(export_window, text=self.t('choose_format'), bg="#1a1a1a",
                fg="white", font=("Arial", 12, "bold")).pack(pady=20)
        
        btn_frame = tk.Frame(export_window, bg="#1a1a1a")
        btn_frame.pack(pady=10)
        
        tk.Button(btn_frame, text=self.t('export_fasta'), command=lambda: [self.export_fasta(), export_window.destroy()],
                 bg="#c62828", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=20, pady=8).pack(side=tk.LEFT, padx=10)
        
        tk.Button(btn_frame, text=self.t('export_pdf'), command=lambda: [self.export_pdf(), export_window.destroy()],
                 bg="#1976d2", fg="white", font=("Arial", 10, "bold"),
                 cursor="hand2", relief=tk.FLAT, bd=0, padx=20, pady=8).pack(side=tk.LEFT, padx=10)
    
    def export_fasta(self):
        """Exporta sequências em formato FASTA"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
            initialfile=f"{self.current_organism.replace(' ', '_')}.fasta"
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    for seq in self.sequences:
                        f.write(f">{seq['id']} {seq['description']}\n")
                        # Quebra sequência em linhas de 60 caracteres
                        sequence = seq['sequence']
                        for i in range(0, len(sequence), 60):
                            f.write(sequence[i:i+60] + '\n')
                        f.write('\n')
                
                messagebox.showinfo('Sucesso', self.t('export_success'))
            except Exception as e:
                messagebox.showerror(self.t('error_title'), f'Erro ao exportar: {str(e)}')
    
    def export_pdf(self):
        """Exporta resultados em formato PDF"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")],
            initialfile=f"{self.current_organism.replace(' ', '_')}_report.pdf"
        )
        
        if filename:
            try:
                doc = SimpleDocTemplate(filename, pagesize=A4, rightMargin=40, leftMargin=40, topMargin=40, bottomMargin=40)
                story = []
                styles = getSampleStyleSheet()
                
                # Título
                title_style = ParagraphStyle(
                    'CustomTitle',
                    parent=styles['Heading1'],
                    fontSize=24,
                    textColor=colors.HexColor('#c62828'),
                    spaceAfter=20,
                    alignment=1  # Center
                )
                
                # Traduz títulos baseado no idioma atual
                lang_titles = {
                    'pt': {'title': 'Relatório NCBI', 'date': 'Data', 'sequences': 'Sequências Encontradas',
                           'sequence': 'Sequência', 'description': 'Descrição', 'size': 'Tamanho',
                           'gc': 'Conteúdo GC', 'composition': 'Composição de Bases', 'alignment': 'Análise de Alinhamento',
                           'comparison': 'Comparação de Sequências', 'type': 'Tipo', 'score': 'Score',
                           'identity': 'Identidade', 'similarity': 'Similaridade', 'gaps': 'Gaps'},
                    'en': {'title': 'NCBI Report', 'date': 'Date', 'sequences': 'Sequences Found',
                           'sequence': 'Sequence', 'description': 'Description', 'size': 'Size',
                           'gc': 'GC Content', 'composition': 'Base Composition', 'alignment': 'Alignment Analysis',
                           'comparison': 'Sequence Comparison', 'type': 'Type', 'score': 'Score',
                           'identity': 'Identity', 'similarity': 'Similarity', 'gaps': 'Gaps'},
                    'es': {'title': 'Informe NCBI', 'date': 'Fecha', 'sequences': 'Secuencias Encontradas',
                           'sequence': 'Secuencia', 'description': 'Descripción', 'size': 'Tamaño',
                           'gc': 'Contenido GC', 'composition': 'Composición de Bases', 'alignment': 'Análisis de Alineamiento',
                           'comparison': 'Comparación de Secuencias', 'type': 'Tipo', 'score': 'Puntuación',
                           'identity': 'Identidad', 'similarity': 'Similitud', 'gaps': 'Huecos'},
                    'fr': {'title': 'Rapport NCBI', 'date': 'Date', 'sequences': 'Séquences Trouvées',
                           'sequence': 'Séquence', 'description': 'Description', 'size': 'Taille',
                           'gc': 'Contenu GC', 'composition': 'Composition des Bases', 'alignment': 'Analyse d\'Alignement',
                           'comparison': 'Comparaison de Séquences', 'type': 'Type', 'score': 'Score',
                           'identity': 'Identité', 'similarity': 'Similarité', 'gaps': 'Écarts'},
                    'de': {'title': 'NCBI-Bericht', 'date': 'Datum', 'sequences': 'Gefundene Sequenzen',
                           'sequence': 'Sequenz', 'description': 'Beschreibung', 'size': 'Größe',
                           'gc': 'GC-Gehalt', 'composition': 'Basenzusammensetzung', 'alignment': 'Alignment-Analyse',
                           'comparison': 'Sequenzvergleich', 'type': 'Typ', 'score': 'Bewertung',
                           'identity': 'Identität', 'similarity': 'Ähnlichkeit', 'gaps': 'Lücken'},
                    'zh': {'title': 'NCBI报告', 'date': '日期', 'sequences': '找到的序列',
                           'sequence': '序列', 'description': '描述', 'size': '大小',
                           'gc': 'GC含量', 'composition': '碱基组成', 'alignment': '比对分析',
                           'comparison': '序列比较', 'type': '类型', 'score': '得分',
                           'identity': '同一性', 'similarity': '相似性', 'gaps': '缺口'},
                    'ru': {'title': 'Отчет NCBI', 'date': 'Дата', 'sequences': 'Найденные последовательности',
                           'sequence': 'Последовательность', 'description': 'Описание', 'size': 'Размер',
                           'gc': 'Содержание GC', 'composition': 'Состав оснований', 'alignment': 'Анализ выравнивания',
                           'comparison': 'Сравнение последовательностей', 'type': 'Тип', 'score': 'Оценка',
                           'identity': 'Идентичность', 'similarity': 'Сходство', 'gaps': 'Пробелы'}
                }
                
                t = lang_titles.get(self.current_language, lang_titles['en'])
                
                story.append(Paragraph(f"{t['title']} - {self.current_organism}", title_style))
                story.append(Spacer(1, 0.1*inch))
                
                # Data e hora
                date_format = '%d/%m/%Y %H:%M' if self.current_language == 'pt' else '%Y-%m-%d %H:%M'
                story.append(Paragraph(f"<b>{t['date']}:</b> {datetime.now().strftime(date_format)}", styles['Normal']))
                story.append(Spacer(1, 0.3*inch))
                
                # Resumo de sequências
                heading2_style = ParagraphStyle('Heading2Custom', parent=styles['Heading2'], textColor=colors.HexColor('#1976d2'))
                story.append(Paragraph(t['sequences'], heading2_style))
                story.append(Spacer(1, 0.15*inch))
                
                for i, seq in enumerate(self.sequences[:10], 1):  # Primeiras 10 sequências
                    story.append(Paragraph(f"<b>{t['sequence']} {i}:</b> {seq['id']}", styles['Normal']))
                    
                    desc_text = seq['description'][:200] + ('...' if len(seq['description']) > 200 else '')
                    story.append(Paragraph(f"<b>{t['description']}:</b> {desc_text}", styles['Normal']))
                    story.append(Paragraph(f"<b>{t['size']}:</b> {seq['length']} bp", styles['Normal']))
                    
                    # Análise de composição
                    sequence = seq['sequence'].upper()
                    counts = Counter(sequence)
                    total = len(sequence)
                    gc = ((counts.get('G', 0) + counts.get('C', 0)) / total * 100) if total > 0 else 0
                    at = ((counts.get('A', 0) + counts.get('T', 0)) / total * 100) if total > 0 else 0
                    
                    story.append(Paragraph(f"<b>{t['gc']}:</b> {gc:.2f}% | AT: {at:.2f}%", styles['Normal']))
                    
                    # Composição detalhada
                    comp_text = f"<b>{t['composition']}:</b> A={counts.get('A',0)} ({counts.get('A',0)/total*100:.1f}%), "
                    comp_text += f"T={counts.get('T',0)} ({counts.get('T',0)/total*100:.1f}%), "
                    comp_text += f"G={counts.get('G',0)} ({counts.get('G',0)/total*100:.1f}%), "
                    comp_text += f"C={counts.get('C',0)} ({counts.get('C',0)/total*100:.1f}%)"
                    story.append(Paragraph(comp_text, styles['Normal']))
                    story.append(Spacer(1, 0.2*inch))
                
                # Adiciona análise de alinhamento se houver dados de comparação
                if hasattr(self, 'last_alignment_data') and self.last_alignment_data:
                    story.append(Paragraph(t['comparison'], heading2_style))
                    story.append(Spacer(1, 0.15*inch))
                    
                    align_data = self.last_alignment_data
                    story.append(Paragraph(f"<b>{t['type']}:</b> {align_data.get('type', 'N/A')}", styles['Normal']))
                    story.append(Paragraph(f"<b>{t['score']}:</b> {align_data.get('score', 'N/A')}", styles['Normal']))
                    story.append(Paragraph(f"<b>{t['identity']}:</b> {align_data.get('identity', 'N/A')}", styles['Normal']))
                    story.append(Paragraph(f"<b>{t['similarity']}:</b> {align_data.get('similarity', 'N/A')}", styles['Normal']))
                    story.append(Paragraph(f"<b>{t['gaps']}:</b> {align_data.get('gaps', 'N/A')}", styles['Normal']))
                    story.append(Spacer(1, 0.2*inch))
                
                doc.build(story)
                messagebox.showinfo(self.t('success_title'), self.t('export_success'))
            except Exception as e:
                messagebox.showerror(self.t('error_title'), f'{self.t("error_export")}: {str(e)}')


def main():
    root = tk.Tk()
    app = NCBISearchApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
