import pydoc
import os
import sys
import shutil


def make_doc(name, parent_dir):
    sys.path.append(parent_dir)
    os.chdir(os.getcwd() + '/docs')
    pydoc.writedoc(name)
    os.chdir(os.path.dirname(os.getcwd()))


def main():
    docs_path = os.getcwd() + '/docs'
    if os.path.exists(docs_path):
        shutil.rmtree(docs_path)

    os.mkdir(docs_path)

    # Docs for tools
    make_doc('tools', os.getcwd() + '/src')

    # Docs for detector
    make_doc('detector', os.getcwd() + '/src/detectors')

    # Docs for scintillator
    make_doc('scintillator', os.getcwd() + '/src/detectors')

    # Docs for shortevent
    make_doc('shortevent', os.getcwd() + '/src/events')

    # Docs for longevent
    make_doc('longevent', os.getcwd() + '/src/events')


if __name__ == '__main__':
    main()
