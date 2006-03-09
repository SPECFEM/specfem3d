#!/usr/bin/env python


from Script import Script
from opal.applications.WebApplication import WebApplication


class CGIScript(Script, WebApplication):

    class Inventory(Script.Inventory, WebApplication.Inventory):
        from pyre.inventory import bool, facility, outputFile
        from opal.components.Login import Login
        actor = facility("actor", factory=Login)
        dump = bool("dump")
        output = outputFile("output")
    
    def main(self, *args, **kwds):
        if self.inventory.dump:
            self.dumpConfiguration()
        self.printWebPage()
        return

    def dumpConfiguration(self):
        configuration = self.retrieveConfiguration()
        print >> self.inventory.output, "\n".join(self.weaver.render(configuration))

    def printWebPage(self):
        import opal.content
        page = opal.content.page()
        head = page.head()
        head.title("SPECFEM3D Basin Solver")
        body = page.body()
        content = body.pageContent()
        main = content.main()
        document = main.document(title='Welcome')
        p = document.paragraph()
        p.text = [
                  '<img src="http://www.gps.caltech.edu/research/jtromp/research/graphics/mpi_slicesThumbnail.gif">',
                  'Hello from the SPECFEM3D Basin Solver!',
                  ]
        if self.inventory.dump:
            p.text.append('<p>Configuration written to %s.' % self.inventory.output.name)
        footer = body.pageFooter()
        self.render(page)


if __name__ == '__main__':
    cgiScript = CGIScript()
    cgiScript.run()


# end of file
