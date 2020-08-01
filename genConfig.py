import polyomino as _mino
import numpy as np
import itertools

def geo_generator(n=1, onlySym=False, exLine=False):
    """
    generate all free geometry configurations
    """
    minos = sorted(_mino.free(_mino.generate(n)), key=_mino.mino_key)
    print('original minos #:', len(minos))
    sorted_minos = []
    if onlySym:
        for mino in minos:
            if mino.symmetry() != '?':
                sorted_minos.append(mino)
        print('symmetric minos #:', len(sorted_minos))
    else:
        sorted_minos = minos


    if exLine:
        exline_minos = []
        for mino in sorted_minos:
            h, w = mino.shape
            mino_list = list(pair for pair in mino)
            # print(mino)

            if h / w == n:
                exline_minos.append(mino)

            isLine = True
            for ind in mino_list:
                num_nhb = 0
                for neighbor in _mino._neighbors(ind):
                    if neighbor in mino_list:
                        num_nhb += 1
                    if num_nhb > 2:
                        isLine = False
                        break

            if not isLine:
                exline_minos.append(mino)
        sorted_minos = exline_minos
        print('exclude redundant line minos #:', len(sorted_minos))

    return sorted_minos

def conn_generator(mino, exSym=True):
    """
    generate all possible connection configs with given geometry config

    @ param: takes in a geometry config
    @ output: all possible connection configs set
    """
    conn_repo = []
    x = [0, 1, 2, 3] # 0 -> rot 0 , 1 -> rot 90 , 2 -> rot 180 , 3 -> rot 270
    ori_conn = [p for p in itertools.product(x, repeat=len(mino))]
    print('...Connection configs num for this geo config:', len(ori_conn))

    if exSym:
        conns = del_symmetric(ori_conn, mino)
    else:
        conns = ori_conn

    for conn in conns:
        conn_config = gen_connlist(mino, conn)
        conn_repo.append(conn_config)
    return conn_repo

def gen_connlist(mino, poses):
    '''
    generate connection list with given geo config with pose info

    @ param: a geometry config
    @ param: a single pose -> list
    @ return: connection list -> list
    '''
    mino_list = list(pair for pair in mino)
    config_dict = {}

    for i in range(len(mino_list)):
        config_dict[mino_list[i]] = (i, poses[i])

    connlist = []
    for cur in mino_list:
        cur_order, cur_pose = config_dict.get(cur)
        config_dict.pop(cur)

        neighbors = _mino._neighbors(cur)
        for neighbor in neighbors:
            if neighbor not in config_dict:
                continue
            ngh_order, ngh_pose = config_dict[neighbor]
            dir12 = tuple(np.asarray(cur) - np.asarray(neighbor))

            relation = side2side(cur_order, ngh_order, dir12, cur_pose, ngh_pose)
            connlist.append(relation)
    # print(connlist)
    return connlist


def del_symmetric(conns, mino):
    result = []
    mapList =[]
    print('sym info',mino.symmetry())
    h, w = mino.shape
    grid_map = np.zeros((h,w))
    sym = mino.symmetry()


    for conn in conns:
        # generate a grid map
        for ind, pose in zip(mino, conn):
            grid_map[ind[0], ind[1]] = pose

        # horizon flip
        if '-' in sym:
            flip_grid = tuple(map(tuple,np.flip(grid_map, 0))) # 0 -> horizon flip  1-> vertical flip
            if flip_grid in mapList:
                continue

        # vertical flip
        if '|' in sym:
            # print('|')
            flip_grid = tuple(map(tuple,np.flip(grid_map, 1))) # 0 -> horizon flip  1-> vertical flip
            if flip_grid in mapList:
                continue

        # Twofold rotational
        if '%' in sym:

            # Fourfold rotational
            if '@' in sym:

                # rot 90 clockwise
                flip_grid = tuple(zip(*grid_map[::-1]))
                if flip_grid in mapList:
                    continue

                # rot 180 clockwise
                flip_grid = tuple(zip(*flip_grid[::-1]))
                if flip_grid in mapList:
                    continue

                # rot 270 clockwise
                flip_grid = tuple(zip(*flip_grid[::-1]))
                if flip_grid in mapList:
                    continue

            # rot 90 clockwise
            flip_grid = tuple(zip(*grid_map[::-1]))

            # rot 180 clockwise
            flip_grid = tuple(zip(*flip_grid[::-1]))

            if flip_grid in mapList:
                continue

        result += [conn]
        mapList += [tuple(map(tuple, grid_map))]
    print(len(conns), len(result))
    return list(result)


def side2side(module1, module2, dir12, pose1, pose2):
    '''
    generate smores connection list by relative position & poses
    '''
    up = {0:'B', 1:'R', 2:'T', 3:'L'}
    down = {0:'T', 1:'L', 2:'B', 3:'R'}
    left = {0:'R', 1:'T', 2:'L', 3:'B'}
    right ={0:'L', 1:'B', 2:'R', 3:'T'}

    if (dir12 == (0, 1)):
        side1 = right[pose1]
        side2 = left[pose2]

    if (dir12 == (0, -1)):
        side1 = left[pose1]
        side2 = right[pose2]

    if (dir12 == (1, 0)):
        side1 = down[pose1]
        side2 = up[pose2]

    if (dir12 == (-1, 0)):
        side1 = up[pose1]
        side2 = down[pose2]

    return [str(module1) + side1, str(module2) + side2]



if __name__ == '__main__':
    minos = geo_generator(n=5, onlySym=True, exLine=True)
    # print('All geo configs num：', len(minos))
    np.random.seed(0)
    random_ind = np.random.randint(0, len(minos))
    random_mino = minos[random_ind]
    print(random_mino)

    all_configs = conn_generator(random_mino)
    print(all_configs)

