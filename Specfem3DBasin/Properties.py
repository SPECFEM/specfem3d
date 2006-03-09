#!/usr/bin/env python


from pyre.inventory.Property import Property


class OutputDir(Property):

    def __init__(self, name, default=None, meta=None, validator=None):
        if default is None:
            default = "."
        Property.__init__(self, name, "directory", default, meta, validator)
        return

    def _cast(self, value):
        from os import makedirs
        from os.path import isdir
        if isinstance(value, basestring):
            if isdir(value):
                pass
            else:
                makedirs(value)
        return value


# end of file
