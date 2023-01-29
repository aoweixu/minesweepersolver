from collections import Counter
import random
import numpy as np
from PIL import Image, ImageGrab
import imagehash
# import matplotlib.pyplot as plt
import pyautogui
from constraint import Problem, ExactSumConstraint

#Setup Parameters:
SQUARE_SIZE = 50 # Cell size in pixels (adjustable)
GREY_RGB_DETECTION_HIGHRES = 198
MINE_TOTAL = 99
BOUNDARY_DETECTION_LIMIT = 0.5 # ratio of boundary to all unopened cells to stop grouping

# Globals
HIGH_RES = False
GRID_WIDTH = 0
GRID_HEIGHT = 0
FIRST_SQUARE_X = 0
FIRST_SQUARE_Y = 0
GUESS_COUNTER = 0

HASH_DICT = {
    0: imagehash.average_hash(Image.open('reference_pictures/empty_80ratio.png')),
    1: imagehash.average_hash(Image.open('reference_pictures/num1_80ratio.png')),
    2: imagehash.average_hash(Image.open('reference_pictures/num2_80ratio.png')),
    3: imagehash.average_hash(Image.open('reference_pictures/num3_80ratio.png')),
    4: imagehash.average_hash(Image.open('reference_pictures/num4_80ratio.png')),
    5: imagehash.average_hash(Image.open('reference_pictures/num5_80ratio.png')),
    6: imagehash.average_hash(Image.open('reference_pictures/num6_80ratio.png')),
    7: imagehash.average_hash(Image.open('reference_pictures/num7_80ratio.png')),
    -1: imagehash.average_hash(Image.open('reference_pictures/mine_80ratio.png')),
    -2: imagehash.average_hash(Image.open('reference_pictures/mine_clicked_80ratio.png')),
    -3: imagehash.average_hash(Image.open('reference_pictures/unopened_80ratio.png')),
    -4: imagehash.average_hash(Image.open('reference_pictures/flag_80ratio.png'))
}
GAME_GRID = None # 2D array of the game grid

def calibrate():
    # Calibrate game grid
    global FIRST_SQUARE_X, FIRST_SQUARE_Y, GAME_GRID, SQUARE_SIZE, HIGH_RES

    # Check screen resolution
    if pyautogui.size()[0] > 1920:
        HIGH_RES = True

    SQUARE_SIZE = 50 if HIGH_RES else 38

    ss = ImageGrab.grab()
    ss = np.array(ss)

    grey_color = [GREY_RGB_DETECTION_HIGHRES,
                  GREY_RGB_DETECTION_HIGHRES,
                  GREY_RGB_DETECTION_HIGHRES]

    Y, X = np.asarray(np.all(ss == grey_color, axis=2)).nonzero()


    Xc, Yc = [], []
    for (xi, yi) in zip(X, Y):
        if (ss[yi, xi-1, :].tolist() != grey_color and
            ss[yi-1, xi, :].tolist() != grey_color and
            ss[yi+SQUARE_SIZE, xi, :].tolist() == grey_color and
            ss[yi+SQUARE_SIZE-1, xi, :].tolist() != grey_color):
            Xc.append(xi)
            Yc.append(yi)

    # Optional additional filtering step depending on screen resolution
    # or game mode (Xc < 300 is easy or intermeidate)
    if HIGH_RES or len(Xc) < 300:
        Xc = Xc[1:]
        Yc = Yc[1:]

    # Manually add last row that is not picked up from initial detection
    row_square_count = int((max(Xc) - min(Xc))/SQUARE_SIZE) + 1
    Ym = max(Yc)

    for i in range(0, row_square_count):
        Xc.append(min(Xc) + i * SQUARE_SIZE)
        Yc.append(Ym + SQUARE_SIZE)

    FIRST_SQUARE_X = Xc[0]
    FIRST_SQUARE_Y = Yc[0]
    row_locs = Counter(Yc).keys()
    col_locs = Counter(Xc).keys()
    GAME_GRID = np.ones(shape=(len(row_locs), len(col_locs))) * -3 # Initialize game grid array of unopened squares

    # For Debug
    # print(FRIST_SQUARE_X)
    # print(FIRST_SQUARE_Y)
    # print(len(Xc))
    # plt.scatter(Xc, Yc, s=2, marker='o')
    # plt.show()
    print('Calibration Complete')

def click_square(x_grid, y_grid):
    # Click a square on the game grid given x and y coordinates in terms of squares
    x_screen, y_screen = convert_coord(x_grid, y_grid)
    print(f'Clicking x={x_screen}({x_grid}) y={y_screen}({y_grid})')
    pyautogui.click(x_screen, y_screen)

def flag_square(x_grid, y_grid):
    # Right click a square on the game grid given x and y coordinates in terms of squares
    x_screen, y_screen = convert_coord(x_grid, y_grid)
    print(f'Flagging x={x_screen}({x_grid}) y={y_screen}({y_grid})')
    pyautogui.rightClick(x_screen, y_screen)

def convert_coord(coord_x, coord_y):
    # Converts grid coordinates to screen (pixel) coordinates
    conv_x = FIRST_SQUARE_X + SQUARE_SIZE * coord_x
    conv_y = FIRST_SQUARE_Y + SQUARE_SIZE * coord_y
    return conv_x, conv_y

def get_square_type(screen, x_grid, y_grid):
    # Try to determine what is in square
    square_capture_ratio = 0.8

    x_screen, y_screen = convert_coord(x_grid, y_grid)
    capture_box = (x_screen,
                   y_screen,
                   x_screen + round(SQUARE_SIZE * square_capture_ratio),
                   y_screen + round(SQUARE_SIZE * square_capture_ratio))

    capture = screen.crop(capture_box)

    # Hash the square image and compare with all square hashes
    hash_square = imagehash.average_hash(capture)

    hash_diff = 100
    square_type = -10
    for square, square_hash in HASH_DICT.items():
        if hash_square - square_hash < hash_diff:
            square_type = square
            hash_diff = hash_square - square_hash

    return square_type

def update_game_grid():
    # Take a screenshot of the game and update the grid array

    global GAME_GRID
    unopened_count = 0
    screen = ImageGrab.grab()

    # Iterate by column then row
    for i in range(0, GAME_GRID.shape[1]):
        for j in range(0, GAME_GRID.shape[0]):
            square_type = get_square_type(screen, i, j)

            # Check for GAME OVER
            if square_type in (-1, -2):
                return 1

            # Count/Check for game solved
            if square_type == -3:
                unopened_count += 1

            GAME_GRID[j][i] = square_type

    if unopened_count:
        return 0 # Game in progress
    else:
        return -1 # Solved

def inside_grid(x_grid, y_grid):
    # Pass grid coordinate through edge detection
    if x_grid < 0 or y_grid < 0 or x_grid >= GAME_GRID.shape[1] or y_grid >= GAME_GRID.shape[0]:
        return -5
    return GAME_GRID[y_grid][x_grid]

def count_flags_around(x_grid, y_grid):
    count = 0
    sq_type_num = -4

    # Count flags in 8 surrounding squares
    if inside_grid(x_grid, y_grid - 1) == sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid) == sq_type_num: count+=1
    if inside_grid(x_grid, y_grid + 1) == sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid) == sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid - 1) == sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid - 1) == sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid + 1) == sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid + 1) == sq_type_num: count+=1

    return count

def count_unopened_around(x_grid, y_grid):
    count = 0
    sq_type_num = (-3, -4)

    # Count unopened squares including flags in 8 surrounding squares
    if inside_grid(x_grid, y_grid - 1) in sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid) in sq_type_num: count+=1
    if inside_grid(x_grid, y_grid + 1) in sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid) in sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid - 1) in sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid - 1) in sq_type_num: count+=1
    if inside_grid(x_grid - 1, y_grid + 1) in sq_type_num: count+=1
    if inside_grid(x_grid + 1, y_grid + 1) in sq_type_num: count+=1

    return count

def is_boundary(x_grid, y_grid):
    # Check whether square is unopened with opened squares in any of surrounding 8 squares

    # Check if square is unopened
    if GAME_GRID[y_grid][x_grid] != -3:
        return False

    boundary = False
    sq_type_num = [1, 2, 3, 4, 5, 6, 7, 8]

    # Edge detection

    # Check if any of 8 neighbour squares are opened
    if inside_grid(x_grid, y_grid - 1) in sq_type_num: boundary = True
    if inside_grid(x_grid - 1, y_grid) in sq_type_num: boundary = True
    if inside_grid(x_grid, y_grid + 1) in sq_type_num: boundary = True
    if inside_grid(x_grid + 1, y_grid) in sq_type_num: boundary = True
    if inside_grid(x_grid - 1, y_grid - 1) in sq_type_num: boundary = True
    if inside_grid(x_grid + 1, y_grid - 1) in sq_type_num: boundary = True
    if inside_grid(x_grid - 1, y_grid + 1) in sq_type_num: boundary = True
    if inside_grid(x_grid + 1, y_grid + 1) in sq_type_num: boundary = True

    return boundary


def random_click(x_grid_1, y_grid_1, x_grid_2, y_grid_2):
    # Click on a random unopened square within the specified rectangle
    # TODO: add handling for infinite looping, i.e. no empty squares found?
    update_game_grid()

    x_grid_rand = random.choice(range(x_grid_1, x_grid_2))
    y_grid_rand = random.choice(range(y_grid_1, y_grid_2))

    while GAME_GRID[y_grid_rand][x_grid_rand] != -3:
        x_grid_rand = random.choice(range(x_grid_1, x_grid_2))
        y_grid_rand = random.choice(range(y_grid_1, y_grid_2))

    print('Guessing...')
    global GUESS_COUNTER
    GUESS_COUNTER += 1
    click_square(x_grid_rand, y_grid_rand)

def attempt_flag_mines():
    # Attempt to flag all guaranteed mines, e.g. square number = surrounding unopened squares

    global GAME_GRID, MINE_TOTAL
    flag_count = 0

    for i in range(0, GAME_GRID.shape[1]):
        for j in range(0, GAME_GRID.shape[0]):

            if GAME_GRID[j][i] > 0:
                surrounding_mines = GAME_GRID[j][i]

                if count_unopened_around(i, j) == surrounding_mines:
                    for ii in range(i-1, i+2):
                        for jj in range(j-1, j+2):
                            # Edge detection
                            if ii < 0 or ii >= GAME_GRID.shape[1] or jj < 0 or jj >= GAME_GRID.shape[0]:
                                continue
                            if GAME_GRID[jj][ii] == -3:
                                flag_square(ii, jj)
                                MINE_TOTAL -= 1
                                GAME_GRID[jj][ii] = -4
    if flag_count:
        print(f'Flagged {flag_count} squares')
    else:
        print('No flaggable squares')

def attempt_open_squares():
    # Attempt to open safe squares guaranteed no mines

    success = False
    for i in range(0, GAME_GRID.shape[1]):
        for j in range(0, GAME_GRID.shape[0]):

            if GAME_GRID[j][i] > 0:
                current_square = GAME_GRID[j][i]
                flags = count_flags_around(i, j)
                unopeneds = count_unopened_around(i, j)

                if current_square == flags and unopeneds > flags:
                    click_square(i, j)
                    success = True

    if success:
        return

    print('No easy cells to open')
    if MINE_TOTAL > 95:
        random_click(0, 0, GAME_GRID.shape[1], GAME_GRID.shape[0])
    else:
        logic_solver()

    # Random guess if no safe way forward
    # print('No safe openable squares')
    # random_click(0, 0, GAME_GRID.shape[1], GAME_GRID.shape[0])

def find_groups(cell_list):
    # Takes a collection of all unopened border cells and attempt to
    # divide them into groups
    # cell_list is a list of (x, y) tuples
    update_game_grid()

    regions = []
    starting_cell = tuple()

    while True:
        # print(f'number of edge cells: {len(cell_list)}') #debug
        if len(cell_list) == 0:
            break

        current_cell_region = []

        # Find cell to process
        starting_cell = random.choice(cell_list)
        # print(f'starting cell: {starting_cell}') #debug
        current_cell_region.append(starting_cell)
        cell_list.remove(starting_cell)

        connections_found = 1

        # Loop all border cells through all current region cells
        while connections_found:
            connections_found = 0
            for cell in cell_list:
                # print(f'Checking cell {cell} against...') #debug
                for region_cell in current_cell_region:
                    # print(f'{region_cell}') #debug
                    if len(cell_list) == 0:
                        break
                    if has_shared_neighbour(cell, region_cell):
                        current_cell_region.append(cell)
                        cell_list.remove(cell)
                        connections_found += 1
                        # print(f'cell removed from list and added to current region. current edge cell count: {len(cell_list)}') #debug
                        break

        regions.append(current_cell_region)

    return regions

def are_adjacent(cell_1, cell_2):
    # Check if two cells are adjacent
    # cell_1 and cell_2 are (x, y) tuples

    x_grid_1 = cell_1[0]
    y_grid_1 = cell_1[1]

    x_grid_2 = cell_2[0]
    y_grid_2 = cell_2[1]

    cell_distance = abs(x_grid_1 - x_grid_2) + abs(y_grid_1 - y_grid_2)

    if cell_1 != cell_2 and cell_distance == 1:
        return True
    else:
        return False

def has_shared_neighbour(cell_1, cell_2):
    # Check if two adjacent cells share an opened neighbour
    # cell_1 and cell_2 are (x, y) tuples

    if not are_adjacent(cell_1, cell_2):
        return False

    x_grid_1 = cell_1[0]
    y_grid_1 = cell_1[1]

    x_grid_2 = cell_2[0]
    y_grid_2 = cell_2[1]

    cell_1_opened_neighbours = []
    for i in range(x_grid_1 - 1, x_grid_1 + 2):
        for j in range(y_grid_1 - 1, y_grid_1 + 2):
            if inside_grid(i, j) > 0:
                cell_1_opened_neighbours.append((i, j))

    for i in range(x_grid_2 - 1, x_grid_2 + 2):
        for j in range(y_grid_2 - 1, y_grid_2 + 2):
            if inside_grid(i, j) > 0 and (i, j) in cell_1_opened_neighbours:
                return True

    return False

def get_regional_opened_cells(cell_region):
    # Finds all opened cells connected to a region of cells

    sq_type_num = [1, 2, 3, 4, 5, 6, 7, 8]
    found_cells = []
    for cell in cell_region:

        for i in range(cell[0] - 1, cell[0] + 2):
            for j in range(cell[1] - 1, cell[1] + 2):
                if (i, j) not in found_cells and inside_grid(i, j) in sq_type_num:
                    found_cells.append((i, j))

    return found_cells

def logic_solver():
    # Attempt to find safe cells using constraint solving
    print('Engaging logic solver..')
    # Find all border pieces
    update_game_grid()
    cell_list = []
    for i in range(GAME_GRID.shape[1]):
        for j in range(GAME_GRID.shape[0]):
            if is_boundary(i, j):
                cell_list.append((i, j))

    # Divide unsolved cells into 

    # Use either one of the two lines below, don't bother finding groups if it doesn't
    # change solving speed significantly, i.e. fast CPU

    # cell_regions = find_groups(cell_list) 
    cell_regions = [cell_list]

    all_solutions = dict()

    for num, region in enumerate(cell_regions):

        problem = Problem()

        # Define variables which are each cell in the group
        sol_dict = dict()
        for cell in region:

            var_name = f'{cell[0]}-{cell[1]}' # X-Y
            problem.addVariable(var_name, [0, 1]) # 0 is safe while 1 is mine

            # Create empty dictionary with variable indexes to store solutions
            sol_dict[var_name] = 0

        # Define constraints, one for each surrounding opened cell
        opened_cells = get_regional_opened_cells(region)
        for cell in opened_cells:

            # Find true mine count
            mine_count = GAME_GRID[cell[1]][cell[0]] - count_flags_around(cell[0], cell[1])

            # Find surrounding unopened cells:
            constraint_vars = []
            for i in range(cell[0] - 1, cell[0] + 2):
                for j in range(cell[1] - 1, cell[1] + 2):
                    if inside_grid(i, j) == -3 and (i, j) in region:
                        constraint_vars.append(f'{i}-{j}')

            problem.addConstraint(ExactSumConstraint(mine_count), constraint_vars)

        # Solve the region and calculate possibilities

        # Add occurances of mine/safe cells in all the solutions
        solutions = problem.getSolutions()
        for solution in solutions:
            for var, value in solution.items():
                sol_dict[var] += value

        # Divide occurances by total number of solutions to get probabilities
        if solutions:
            for var, value in sol_dict.items():
                sol_dict[var] = value / len(solutions)

        # Debug
        # print(f'Region {num} solutions:')
        # print(dict(sorted(sol_dict.items(), key=lambda x: x[1])))

        # Merge solutions into one dictionary
        all_solutions = all_solutions | sol_dict

    print(f'All region solutions: {all_solutions}') #debug
    if all_solutions:
        process_solutions(all_solutions)
    return

def process_solutions(solutions_dict):
    global GUESS_COUNTER

    # First try to extract all 100% safe or mine cells
    safe_cells = []
    mine_cells = []

    for cell_key, mine_prob in solutions_dict.items():

        cell_x = int(cell_key.split('-')[0])
        cell_y = int(cell_key.split('-')[1])

        if mine_prob == 1:
            mine_cells.append((cell_x, cell_y))
        if mine_prob == 0:
            safe_cells.append((cell_x, cell_y))

    # Call click and flag methods on the two lists
    print(f'Logic solving complete with {len(safe_cells) + len(mine_cells)} safe/mine cells found')
    for cell in mine_cells:
        flag_square(cell[0], cell[1])
    for cell in safe_cells:
        click_square(cell[0], cell[1])

    if safe_cells or mine_cells:
        return

    # If no guaranteed step then guess the highest probability move
    highest_prob_mine = 0
    mine_cell = tuple()
    highest_prob_safe = 0
    safe_cell = tuple()

    for cell_key, mine_prob in solutions_dict.items():
        cell_x = int(cell_key.split('-')[0])
        cell_y = int(cell_key.split('-')[1])
        if mine_prob > highest_prob_mine:
            highest_prob_mine = mine_prob
            mine_cell = ((cell_x, cell_y))
        if (1 - mine_prob) > highest_prob_safe:
            highest_prob_safe = 1 - mine_prob
            safe_cell = ((cell_x, cell_y))

    if highest_prob_mine > highest_prob_safe:
        print(f'Guessing mine with probability of {round(highest_prob_mine * 100)}%')
        flag_square(mine_cell[0], mine_cell[1])
    else:
        print(f'Guessing safe with probability of {round(highest_prob_safe * 100)}%')
        click_square(safe_cell[0], safe_cell[1])

    GUESS_COUNTER += 1
    return

def main():
    #Start by game board calibration
    calibrate()
    print(GAME_GRID.shape)
    input('Press enter to continue:')

    # Debug // Print grids

    # boundary_list = []
    # boundary_grid = np.zeros(shape=GAME_GRID.shape)
    # for i in range(0, GAME_GRID.shape[1]):
    #     for j in range(0, GAME_GRID.shape[0]):
    #         if is_boundary(i, j):
    #             boundary_grid[j][i] = 1
    #             boundary_list.append((i, j))
    #         else:
    #             boundary_grid[j][i] = 0


    # edge_groups = find_groups()
    # for idx, group in enumerate(edge_groups):
    #     for cell in group:
    #         GAME_GRID[cell[1]][cell[0]] = 99 - (idx * 10)
    # np.savetxt('game_grid.csv', GAME_GRID, delimiter=',')
    # np.savetxt('edge_test.csv', boundary_grid, delimiter=',')

    # START GAME

    random_click(0, 0, GAME_GRID.shape[1], GAME_GRID.shape[0])

    game_over = 0
    step_counter = 0
    while not game_over:

        step_counter+=1
        attempt_flag_mines()
        attempt_open_squares()
        game_over = update_game_grid()

        # Debug

        # np.savetxt(f'grids\step{step_counter}_grid.csv', GAME_GRID, delimiter=',')

    result = 'Loss' if game_over == 1 else 'Win'
    print(f'Game Over. Result: {result}')
    print(f'Steps: {step_counter}')
    print(f'Guesses: {GUESS_COUNTER}')

if __name__ == "__main__":
    main()
    