# fmt: off

def pytest_generate_tests(metafunc):
    if 'seed' in metafunc.fixturenames:
        seeds = metafunc.config.getoption('seed')
        if len(seeds) == 0:
            seeds = [0]
        else:
            seeds = list(map(int, seeds))
        metafunc.parametrize('seed', seeds)


def pytest_addoption(parser):
    parser.addoption(
        '--seed',
        action='append',
        default=[],
        help='add a seed for tests where random number generators'
        ' are involved. This option can be applied more'
        ' than once.',
    )
