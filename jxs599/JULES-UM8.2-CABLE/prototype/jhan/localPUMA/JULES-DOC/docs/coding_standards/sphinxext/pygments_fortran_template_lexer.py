# -*- coding: utf-8 -*-
"""
    Pygments lexer for the Fortran templates in the JULES Coding Standards
"""

import re

from pygments.lexer import RegexLexer, using
from pygments.token import *
from pygments.lexers.compiled import FortranLexer


class FortranTemplateLexer(FortranLexer):
    """Pygments lexer for the Fortran templates used in the JULES Coding Standards"""
    name = 'Fortran Template'
    aliases = ['fortran_template']
    
    flags = re.IGNORECASE | re.DOTALL

    tokens = {
        'root' : [
            # We deal with comments starting at the start of a line separately
            (r'!', Comment.Single, 'comment'),
            
            # Anything up to an angle bracket or a new line uses the Fortran lexer
            # We want to include the new line in the lexed string if it is there
            (r'[^<\n]+(\n)?', using(FortranLexer)),
            
            # Anything between angle brackets is special
            (r'<[^>]+>', Comment.Special),
        ],
        
        'comment' : [
            # Anything up to an angle bracket or a new line is comment text
            (r'[^<\n]+', Comment.Single),
            
            # Anything between angle brackets is special
            (r'<[^>]+>', Comment.Special),
            
            # A new line takes us out of comment mode
            (r'\n', Text, '#pop'),
        ]
    }


def setup(app):
    """Makes the lexer known to Pygments through Sphinx"""
    app.add_lexer('fortran_template', FortranTemplateLexer())