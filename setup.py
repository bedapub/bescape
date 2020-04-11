from distutils.core import setup
setup(
  name = 'bescape',         # How you named your package folder (MyLib)
  packages = ['bescape'],   # Chose the same as "name"
  version = '0.2',      # Start with a small number and increase it with every change you make
  license='gpl-3.0',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'BESCA proportion estimator - BESCAPE is a cell deconvolution package. The user can specify a custom basis vector, as well as the preferred deconvolution method. Thus it allows us to detach the deconvolution algorithm from the underlying basis vector it originally comes packaged with.',   # Give a short description about your library
  author = 'Miro Phan',                   # Type in your name
  author_email = 'phanmir@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/bedapub/bescape',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/bedapub/bescape/releases/tag/v0.2',    # I explain this later on
  keywords = ['bescape', 'deconvolution'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'spython',
          'docker',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.7'
  ],
)
