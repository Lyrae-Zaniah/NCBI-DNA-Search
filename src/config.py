"""
Configurações globais do projeto
"""

# Configuração de email NCBI (será substituído pela interface)
NCBI_EMAIL = ""

# Configurações de SSL
import ssl
import certifi
import urllib.request

try:
    ssl_context = ssl.create_default_context(cafile=certifi.where())
    https_handler = urllib.request.HTTPSHandler(context=ssl_context)
    opener = urllib.request.build_opener(https_handler)
    urllib.request.install_opener(opener)
except:
    ssl._create_default_https_context = ssl._create_unverified_context

# Configurações de alinhamento
ALIGNMENT_CONFIG = {
    'match_score': 2,
    'mismatch_score': -1,
    'gap_open': -2.0,
    'gap_extend': -2.0
}

# Configurações de UI
UI_COLORS = {
    'bg_dark': '#1a1a1a',
    'bg_medium': '#2d2d2d',
    'bg_light': '#3d3d3d',
    'accent_red': '#c62828',
    'accent_blue': '#1976d2',
    'text_white': 'white',
    'text_gray': '#999',
    'border': '#555'
}

# Configurações de fonte
FONTS = {
    'title': ('Arial', 24, 'bold'),
    'heading': ('Arial', 14, 'bold'),
    'normal': ('Arial', 10),
    'code': ('Courier New', 9)
}
