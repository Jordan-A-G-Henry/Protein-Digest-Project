# -*- coding: utf-8 -*-
"""Main Controller"""

from tg import expose, flash, require, url, lurl
from tg import request, redirect, tmpl_context
from tg.i18n import ugettext as _, lazy_ugettext as l_
from tg.exceptions import HTTPFound

from spliceapp.lib.base import BaseController
from spliceapp.controllers.error import ErrorController

__all__ = ['RootController']

from Splice_Mod import *
from spliceapp.model.Splice_App import *

import tw2.forms as twf

#[1]

class SearchForm(twf.Form):
    class child(twf.TableLayout):
        enzyme = twf.SingleSelectField(label="Enzyme", options = ['Trypsin','Lys_C','Lys_N','CNBr','Arg_C','Asp_C'], prompt_text=None)
        rawseq = twf.TextField(label="Raw Sequence")
        min_l = twf.TextField(label="Min. Splice Legnth")
        max_l = twf.TextField(label="Max. Splice Legnth")
        min_w = twf.TextField(label="Min. Splice Weight")
        max_w = twf.TextField(label="Max. Splice Weight")
        mc = twf.TextField(label="Missed Cleavages")

        css_class = 'table'
        attrs = {'style': 'width: 600px;'}


    action = '/dice'
    submit = twf.SubmitButton(value="Search")



class RootController(BaseController):
    @expose('spliceapp.templates.index')
    def index(self):
        """Handle the front-page."""
        return dict(title="Splice-App!", form=SearchForm)

#[2]
    @expose(template="spliceapp.templates.dice")
    def dice(self,enzyme,rawseq, min_l, max_l, min_w, max_w, mc): 
        if min_l == '':
            min_l = 0
        if max_l == '':
            max_l = 1000
        if min_w == '':
            min_w = 0
        if max_w == '':
            max_w = 100000
        if mc == '':
            mc = 0
        
        """Return the table of splice results."""
        Query = Display(enzyme,rawseq, int(min_l), int(max_l), int(min_w), int(max_w), int(mc))
        return dict(title="Splice Results", params=Query.params, splice=Query.table)


