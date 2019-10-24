from distutils.core import setup

setup(name='pygeos',
      version='0.4.0',
      description='GEOS front-end preprocessing package',
      author='Chris Sherman',
      author_email='sherman27@llnl.gov',
      packages=['pygeos', 'pygeos.tests'],
      entry_points={'console_scripts': ['pygeos = pygeos.__main__:main',
                                        'test_pygeos = pygeos.tests.test_manager:run_unit_tests']})

