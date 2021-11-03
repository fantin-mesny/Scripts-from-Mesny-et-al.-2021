import argparse
import sys
from lxml import etree

def get_params(argv):
    parser = argparse.ArgumentParser(description='retrieve .sh file to download from JGI')
    parser.add_argument('-xml', '--xml', help="XML file directory", required=True)
    parser.add_argument('-o', '--o', help="output directory", required=True)
    a = parser.parse_args()
    return a


if __name__ == '__main__':
    a = get_params(sys.argv[1:])
    toDownload=['Files/Assembly/Assembled scaffolds (masked)','Files/Annotation/Filtered Models ("best")/Genes']
    root = etree.fromstring(open(a.xml).read())
    tree = etree.ElementTree(root)
    with open(a.o+'/download.sh', 'w+') as outp:
        outp.write('#!/bin/bash \n')
        for ele in toDownload:
            s=''
            for fold in ele.split('/'):
                s+="/folder[@name='"+fold+"']"
            subtree=tree.xpath("/organismDownloads"+s)[0]

            for e in subtree.iter():
                nodes = root.xpath(tree.getpath(e))
                for node in nodes:
                    if list(node.attrib.keys())!=['name']:
                        outp.write('curl https://genome.jgi.doe.gov%s -b cookies > %s/%s\n' % (node.attrib['url'], a.o, node.attrib['filename']))

