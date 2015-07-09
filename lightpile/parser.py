# -*- coding: utf-8 -*-
"""
Parse the user-generated assignment (*lmc) file.
"""
from __future__ import unicode_literals
from __future__ import print_function
import pyparsing as pp
import re
import sys

import lightpile

__author__ = "Richard Pfeifer"


class Parser(object):
    """
    The Parser builds a dictionary from the assignment file. It identifies and
    replaces all numbers, keywords and identifiers. It does not build any
    objects. It checks filenames for plausability but not for file existence.
    It checks all sections for lexic and syntax but not for validity and
    completeness of the tasks (i.e. it does not check if all materials in the
    stack section are defined in the material section).

    pyparsing is a parser for contextfree grammar. Introducing context is done
    by multi-step parsing. First we parse the input without context and
    identify sections of same context. Then parse the sections with their
    own grammar, each being context free.

    -----GRAMMAR-sections-------
    <grammar> ::= <section> | <section> <grammar>
    <section> ::= <sectionhead> <sectionbody>
    <sectionhead> ::= [ SectionName ] LineEnd
    <sectionbody> ::= <sectionline> | <sectionline> <sectionbody>
    <sectionline> ::= ~<sectionbody> ~StringEnd RestOfLine LineEnd

    -----GRAMMAR material section-------
    SectionName: materials
    <sectionline> ::= <material_line>

    <material_line> ::= t_word t_real |
                        t_word t_real t_real |
                        t_word t_enclosed_string

    -----GRAMMAR stack section------------
    SectionName: stack
    <sectionline> ::= <planarlayer_line> |
                      <emitter_line>

"""
    def __init__(self):

        # codec used for encoding of usermessages
        self.codec = lightpile.codec

        EOL = pp.LineEnd()
        SOL = pp.LineStart().leaveWhitespace()
        blankline = SOL + EOL
        comment = '#' + pp.restOfLine + EOL

        self.comment_stripper = comment.setParseAction(pp.replaceWith("\n"))
        self.blankline_stripper = blankline.setParseAction(pp.replaceWith(""))

        # filegrammar
        ws = ' \t'
        standard_chars = pp.printables.replace('#', '')
        text = pp.OneOrMore(pp.White(ws) |
                            pp.quotedString |
                            pp.Word(standard_chars))
        text.setParseAction(lambda tokens: ''.join(tokens))

        lbracket = pp.Literal("[").suppress()
        rbracket = pp.Literal("]").suppress()
        # identifiers have to start with a unicode-letter and can continue
        # with any number of unicode-letters or any of
        # ':', '%', '+', '-', '_', '.', ','
        self.identifier = pp.Regex(r'[^\W\d]([^\W]|[%:_,\+\-\.])*', re.U)
        sectionhead = (lbracket + self.identifier + rbracket +
                       EOL.suppress())
        sectionline = ~sectionhead + ~pp.StringEnd() + pp.restOfLine + EOL
        sectionlines = pp.ZeroOrMore(sectionline)
        sectionlines.setParseAction(lambda tokens: ''.join(tokens))

        self.filegrammar = pp.dictOf(sectionhead, sectionlines)

        self._init_sectiongrammars()

    def _init_sectiongrammars(self):

        EOL = pp.LineEnd()
        floatNumber = pp.Regex(r'-?\d+(\.\d*)?([eE]-?\d+)?').setParseAction(
            lambda f: float(f[0]))
        intNumber = pp.Regex(r'-?[0-9]+').setParseAction(lambda i: int(i[0]))
        number = (floatNumber | intNumber)

        # kw_spectralrange = pp.Keyword("spectralrange")
        kw_spectralpoint = pp.Keyword("spectralpoint")
        kw_angularrange = pp.Keyword("angularrange")

        spectralquantities = "wavelength frequency energy"
        spectralunits = "nm"
        angularquantities = "kt u"
        angularunits = "1/nm None"
        scalespecifier = "linear log"

        spectralpoint_line = pp.Group(
            kw_spectralpoint.suppress() +
            pp.oneOf(spectralquantities).setResultsName("quantity") +
            pp.oneOf(spectralunits).setResultsName("unit") +
            number.setResultsName("value") +
            EOL.suppress()
        ).setResultsName("spectralpoint")

        sampling_option = (intNumber | 'adaptive')
        angularrange_line = pp.Group(
            kw_angularrange.suppress() +
            pp.oneOf(angularquantities).setResultsName("quantity") +
            pp.oneOf(angularunits).setResultsName("unit") +
            number.setResultsName("begin") +
            number.setResultsName("end") +
            sampling_option.setResultsName("sampling_points") +
            pp.Optional(
                pp.oneOf(scalespecifier), default=None
            ).setResultsName("scale") +
            EOL.suppress()
        ).setResultsName("angularrange")

        # section: materials
        materialname = self.identifier
        opticalconstant = number + pp.Optional(number)
        opticalconstant = opticalconstant.setResultsName("nk")

        def nkToComplex(toks):
            if len(toks) == 2:
                return complex(toks[0], toks[1])
            if len(toks) == 1:
                return complex(toks[0])
            else:
                raise ValueError()
        opticalconstant.setParseAction(nkToComplex)
        dispersivematerial_line = pp.Group(
            materialname.setResultsName("name") +
            opticalconstant +
            EOL.suppress()
        )
        nondispersivematerial_line = pp.Group(
            materialname.setResultsName("name") +
            pp.QuotedString('"').setResultsName("nk") +
            EOL.suppress()
        )
        material_line = dispersivematerial_line | nondispersivematerial_line
        material_sectionbody = pp.OneOrMore(material_line)

        # section: stack
        planar = pp.Group(materialname("planarmaterial") + number("thickness"))
        planar_line = planar + EOL.suppress()
        emittername = self.identifier("emittername")
        emitterposition_line = pp.Group(emittername) + EOL.suppress()
        stack_line = (planar_line | emitterposition_line)
        thinfilm_line = stack_line + pp.FollowedBy(stack_line)
        stack_sectionbody = (planar_line +
                             pp.ZeroOrMore(thinfilm_line) +
                             planar_line)

        # section: emitters
        emittername = self.identifier.setResultsName("name")
        # pl_file = pp.QuotedString('"').setResultsName("pl_file")
        orientation = pp.oneOf("vert hor iso").setResultsName("orientation")
        # emitterweight = number.setResultsName("emitterweight")
        iqe = number.setResultsName("iqe")
        emitter_line = pp.Group(
            emittername + orientation + iqe + EOL.suppress()
        )
        emitters_sectionbody = pp.OneOrMore(emitter_line)

        # section: dipolestudy
        dipolestudy_sectionbody = pp.Each([spectralpoint_line,
                                           angularrange_line])

        # section: LOE
#        outcouplinglimits_line = pp.Group(
#                kw_outcoupling.suppress() +
#                pp.oneOf("top bot").setResultsName("halfspacedirection") +
#                number.setResultsName("max") + pp.Optional(number) +
#                EOL.suppress()
#                )
#        outcouplinglimits_line = \
#                    outcouplinglimits_line.setResultsName('outcouplinglimits')
#        angularmax_line = pp.Group(
#                                kw_umax.suppress() + number + EOL.suppress()
#                                   )
#        angularmax_line = angularmax_line.setResultsName('angularmax')

#        LOE_sectionbody = pp.Each([spectralrange_line,
#                                   outcouplinglimits_line,
#                                   angularmax_line])

        # section: graph_x_a
        xlimits_line = (
            pp.Keyword("xlim").suppress() +
            (number | pp.Keyword("auto")).setResultsName("xlim_left") +
            (number | pp.Keyword("auto")).setResultsName("xlim_right") +
            EOL.suppress()
        )
        ylimits_line = (
            pp.Keyword("ylim").suppress() +
            (number | pp.Keyword("auto")).setResultsName("ylim_bottom") +
            (number | pp.Keyword("auto")).setResultsName("ylim_top") +
            EOL.suppress()
        )
        xscale_line = (
            pp.Keyword("xscale").suppress() +
            (
                pp.Keyword("linear") | pp.Keyword("log")
            ).setResultsName("xscale") +
            EOL.suppress()
        )
        yscale_line = (
            pp.Keyword("yscale").suppress() +
            (
                pp.Keyword("linear") | pp.Keyword("log")
            ).setResultsName("yscale") +
            EOL.suppress()
        )
        lines_line = (
            pp.Keyword("lines").suppress() +
            pp.oneOf("all total").setResultsName("lines") +
            EOL.suppress()
        )
        externaldata_line = pp.Group(
            pp.Keyword("externaldata").suppress() +
            pp.QuotedString('"').setResultsName("path") +
            pp.Optional(
                pp.QuotedString('"'),
                default=None).setResultsName("label")
        ).setResultsName("externaldata")

        graph_option_line = (xlimits_line | ylimits_line | lines_line |
                             xscale_line | yscale_line | externaldata_line)
        graph_x_a_sectionbody = pp.OneOrMore(graph_option_line)

        self.sectiongrammars = {'materials': material_sectionbody,
                                'emitters': emitters_sectionbody,
                                'stack': stack_sectionbody,
                                'dipolestudy': dipolestudy_sectionbody,
                                # 'LOE': LOE_sectionbody,
                                'graph_f_a': graph_x_a_sectionbody,
                                'graph_p_a': graph_x_a_sectionbody
                                }

    def _as_dict(self, parse_results):
        """
        Recursively converts the ParseResults object to a dict.
        """
        d = parse_results.asDict()
        for i in d:
            if type(d[i]) == pp.ParseResults:
                d[i] = self._as_dict(d[i])
        return d

    def _as_list_of_dicts(self, parse_result):
        """
        Converts the ParseResult object to a list of dicts.
        """
        return [self._as_dict(res) for res in parse_result]

    def parse(self, usertext_string):
        """
        Convert (unicode) usertext to dictionary acceptable for 'taskcreator'.

        *usertext_string* is stripped of comments and blank lines. After
        parsing a dictionary is returned. It consists of sectionnames as keys
        and nested dicts or a list of nested dicts as values.

        Args:
            usertext_string Unicode     The userinput (from a user generated
                                        lmc-file)
        Exceptions:
            TypeError       if bytestring is given instead of unicode string.
        """
        if not isinstance(usertext_string, unicode):
            msg = "Parser.parse expects a unicode string."
            raise TypeError(msg)

        uncommented_text = \
            self.comment_stripper.transformString(usertext_string)
        cleaned_text = \
            self.blankline_stripper.transformString(uncommented_text)
        parsed_sections = self.filegrammar.parseString(
            cleaned_text,
            parseAll=True
        )
        sectiondict = {}
        for sectionname, sectionbody in parsed_sections.asDict().iteritems():
            if sectionname not in self.sectiongrammars.keys():
                raise pp.ParseFatalException(
                    "Invalid sectionname '{0}'".format(sectionname))
            grammar = self.sectiongrammars[sectionname]
            try:
                r = grammar.parseString(sectionbody, parseAll=True)
            except pp.ParseException as e:
                msg = "Error parsing section [{0}]\n".format(sectionname)
                sectionlines_list = sectionbody.split('\n')
                if sectionlines_list[-1] == '':
                    # remove last (empty) element
                    sectionlines_list = sectionlines_list[:-1]
                msg = []
                for nr, line in enumerate(sectionlines_list):
                    msg.append("{0}\t{1}\n".format(nr + 1, line))
                msg.append(e.__repr__())
                msg.append("\nPossible cause: Unknown or mistyped keyword ")
                msg.append("after given expected end of text.")
                msg = ''.join(msg)
                print(msg.encode(lightpile.codec), file=sys.stderr)
                sys.exit(1)
            else:
                if sectionname in ["materials", "stack", "emitters"]:
                    r = self._as_list_of_dicts(r)
                elif sectionname in ["dipolestudy", "graph_f_a", "graph_p_a"]:
                    r = self._as_dict(r)
                sectiondict.update({sectionname: r})
        return sectiondict
