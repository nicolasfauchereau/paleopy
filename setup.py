from distutils.core import setup

setup(name='paleopy',
      version='0.1.0',
      description="""
      classes and methods to calculates and plot analog composites from
      a paleoclimate proxy or an ensemble of paleoclimate proxies
      """,
      author='Nicolas Fauchereau',
      author_email='nicolas.fauchereau@gmail.com',
      url='https://github.com/nicolasfauchereau/paleopy',
      packages=['paleopy',
                'paleopy.core',
                'paleopy.plotting',
                'paleopy.markov',
                'paleopy.utils',
                ],
      license='MIT',
      platforms='any',
      classifiers=[
          'License :: MIT',
          'Programming Language :: Python :: 3',
      ],
      long_description="""

A library implementing Python classes, methods and functions
for paleoclimate reconstructions from a proxy or ensemble of proxy
by the methof of analogs.

Reference
=============
Lorrey, A.M., Fauchereau, N., Stanton, C., Chappell, P.R., Phipps, S.J., Mackintosh, A., Renwick, J.A., and Fowler, A.M. 2013.
The Little Ice Age climate of New Zealand reconstructed from Southern Alps cirque glaciers: a synoptic type approach.
Climate Dynamics, June 2014, Volume 42, Issue 11-12, pp 3039-3060.

Contact
=============

If you have any questions or comments about paleopy, please feel free to contact me via
e-mail: nicolas.fauchereau@gmail.com
or Twitter: https://twitter.com/nfauchereau

This project is hosted at https://github.com/nicolasfauchereau/paleopy

""",
    )
