#Copyright 2011 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
The todomod sphinx extension adds the `todomodule` directive which will show all
the ``todo`` directives in a given module all in one place in the sphinx
documentation.  Like the :mod:`sphinx.ext.todo` extension, it will not be shown 
for non-release versions.
"""

from sphinx.ext.todo import Todo,todo_node,nodes
from sphinx.pycode import ModuleAnalyzer,PycodeError


class TodoModule(Todo):
    required_arguments = 1
    has_content = True
    
    def run(self):        
        try:
            modfn = ModuleAnalyzer.for_module(self.arguments[0]).srcname
        except PycodeError,e:
            warnstr = "can't find module %s for todomodule: %s"%(self.arguments[0],e)
            return [self.state.document.reporter.warning(warnstr,lineno=self.lineno)]
        
        todolines = []
        with open(modfn) as f:
            for l in f:
                if l.startswith('#TODO'):
                    todolines.append(l)
                    
        todoreses = []
        for tl in todolines:
            text = tl.replace('#TODO:','').replace('#TODO','').strip()
            env = self.state.document.settings.env

            targetid = "todo-%s" % env.index_num
            env.index_num += 1
            targetnode = nodes.target('', '', ids=[targetid])
            
            td_node = todo_node(text)
            
            title_text = _('Module Todo')
            textnodes, messages = self.state.inline_text(title_text, self.lineno)
            td_node += nodes.title(title_text, '', *textnodes)
            td_node += messages
            if 'class' in self.options:
                classes = self.options['class']
            else:
                classes = ['admonition-' + nodes.make_id(title_text)]
            td_node['classes'] += classes
            
            td_node.append(nodes.paragraph(text,text))
            td_node.line = self.lineno
            
            todoreses.append(targetnode)
            todoreses.append(td_node)
        return todoreses
    

def setup(app):
    app.add_directive('todomodule', TodoModule) #add this directive to document TODO comments in the root of the module
