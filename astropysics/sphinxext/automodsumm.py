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
The automodsum sphinx extension ...
"""

from sphinx.ext.autosummary import Autosummary

class AutoModSumm(Autosummary):
    #inherit from AutoSummary
    def run(self):
        print 'in run'
        return Autosummary.run(self)

def setup(app):
    app.setup_extension('sphinx.ext.autosummary')
    app.add_directive('autosummary', AutoModSumm)