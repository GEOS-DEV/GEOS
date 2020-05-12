from distutils.core import setup

setup(name='pygeos',
      version='0.4.2',
      description='GEOS front-end preprocessing package',
      author='Chris Sherman',
      author_email='sherman27@llnl.gov',
      packages=['pygeos', 'pygeos.tests'],
      entry_points={'console_scripts': ['pygeos = pygeos.__main__:main',
                                        'format_xml = pygeos.__main__:format_xml',
                                        'test_pygeos = pygeos.tests.test_manager:run_unit_tests']},
      install_requires=['lxml>=4.5.0', 'numpy'])

