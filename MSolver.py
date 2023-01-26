import numpy as np
from PIL import Image, ImageGrab
import imagehash
import matplotlib.pyplot as plt
import pyautogui
import random
from collections import Counter

#Setup Parameters:
SQUARE_SIZE = 50
GREY_RGB_DETECTION_HIGHRES = 198

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

    # Optional additional filtering step depending on screen resolution or game mode (Xc < 300 is easy or intermeidate)
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

def count_flags_around(x_grid, y_grid):
    count = 0
    sq_type_num = -4
    # Edge detection
    top = True if y_grid == 0 else False
    left = True if x_grid == 0 else False
    bottom = True if y_grid == GAME_GRID.shape[0] - 1 else False
    right = True if x_grid == GAME_GRID.shape[1] - 1 else False

    # Count flags in 8 surrounding squares
    if not top and GAME_GRID[y_grid - 1][x_grid] == sq_type_num: count+=1
    if not left and  GAME_GRID[y_grid][x_grid - 1] == sq_type_num: count+=1
    if not bottom and GAME_GRID[y_grid + 1][x_grid] == sq_type_num: count+=1
    if not right and GAME_GRID[y_grid][x_grid + 1] == sq_type_num: count+=1
    if not top and not left and GAME_GRID[y_grid - 1][x_grid - 1] == sq_type_num: count+=1
    if not top and not right and GAME_GRID[y_grid - 1][x_grid + 1] == sq_type_num: count+=1
    if not bottom and not left and GAME_GRID[y_grid + 1][x_grid - 1] == sq_type_num: count+=1
    if not bottom and not right and GAME_GRID[y_grid + 1][x_grid + 1] == sq_type_num: count+=1

    return count

def count_unopened_around(x_grid, y_grid):
    count = 0
    sq_type_num = -2
    # Edge detection
    top = True if y_grid == 0 else False
    left = True if x_grid == 0 else False
    bottom = True if y_grid == GAME_GRID.shape[0] - 1 else False
    right = True if x_grid == GAME_GRID.shape[1] - 1 else False

    # Count unopened squares including flags in 8 surrounding squares
    if not top and GAME_GRID[y_grid - 1][x_grid] < sq_type_num: count+=1
    if not left and  GAME_GRID[y_grid][x_grid - 1] < sq_type_num: count+=1
    if not bottom and GAME_GRID[y_grid + 1][x_grid] < sq_type_num: count+=1
    if not right and GAME_GRID[y_grid][x_grid + 1] < sq_type_num: count+=1
    if not top and not left and GAME_GRID[y_grid - 1][x_grid - 1] < sq_type_num: count+=1
    if not top and not right and GAME_GRID[y_grid - 1][x_grid + 1] < sq_type_num: count+=1
    if not bottom and not left and GAME_GRID[y_grid + 1][x_grid - 1] < sq_type_num: count+=1
    if not bottom and not right and GAME_GRID[y_grid + 1][x_grid + 1] < sq_type_num: count+=1

    return count

def is_boundary(x_grid, y_grid):
    # Check whether square is unopened with opened squares in any of surrounding 8 squares

    # Check if square is unopened
    if GAME_GRID[y_grid][x_grid] != -3:
        return False

    boundary = False
    sq_type_num = [1, 2, 3, 4, 5, 6, 7, 8]

    # Edge detection
    top = True if y_grid == 0 else False
    left = True if x_grid == 0 else False
    bottom = True if y_grid == GAME_GRID.shape[0] - 1 else False
    right = True if x_grid == GAME_GRID.shape[1] - 1 else False

    # Check if any of 8 neighbour squares are opened
    if not top and GAME_GRID[y_grid - 1][x_grid] in sq_type_num: boundary = True
    if not left and  GAME_GRID[y_grid][x_grid - 1] in sq_type_num: boundary = True
    if not bottom and GAME_GRID[y_grid + 1][x_grid] in sq_type_num: boundary = True
    if not right and GAME_GRID[y_grid][x_grid + 1] in sq_type_num: boundary = True
    if not top and not left and GAME_GRID[y_grid - 1][x_grid - 1] in sq_type_num: boundary = True
    if not top and not right and GAME_GRID[y_grid - 1][x_grid + 1] in sq_type_num: boundary = True
    if not bottom and not left and GAME_GRID[y_grid + 1][x_grid - 1] in sq_type_num: boundary = True
    if not bottom and not right and GAME_GRID[y_grid + 1][x_grid + 1] in sq_type_num: boundary = True

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

    global GAME_GRID
    flag_count = 0

    for i in range(0, GAME_GRID.shape[1]):
        for j in range(0, GAME_GRID.shape[0]):
            print(f'square i:{i} j:{j}') # Debug
            if GAME_GRID[j][i] > 0:
                surrounding_mines = GAME_GRID[j][i]
                print(f'unopened={count_unopened_around(i,j)}') #debug
                if count_unopened_around(i, j) == surrounding_mines:
                    for ii in range(i-1, i+2):
                        for jj in range(j-1, j+2):
                            # Edge detection
                            if ii < 0 or ii >= GAME_GRID.shape[1] or jj < 0 or jj >= GAME_GRID.shape[0]:
                                continue
                            if GAME_GRID[jj][ii] == -3:
                                flag_square(ii, jj)
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

    # Random guess if no safe way forward
    print('No safe openable squares')
    random_click(0, 0, GAME_GRID.shape[1], GAME_GRID.shape[0])

def main():
    #Start by game board calibration
    calibrate()
    print(GAME_GRID.shape)
    input('Press enter to continue:')

    update_game_grid()

    # Debug // Print grids

    # np.savetxt('game_grid.csv', GAME_GRID, delimiter=',')

    # boundary_grid = np.zeros(shape=GAME_GRID.shape)
    # for i in range(0, GAME_GRID.shape[1]):
    #     for j in range(0, GAME_GRID.shape[0]):
    #         if is_boundary(i, j):
    #             boundary_grid[j][i] = 1
    #         else:
    #             boundary_grid[j][i] = 0

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

        np.savetxt(f'grids\step{step_counter}_grid.csv', GAME_GRID, delimiter=',')

    result = 'Loss' if game_over == 1 else 'Win'
    print(f'Game Over. Result: {result}')
    print(f'Steps: {step_counter}')
    print(f'Guesses: {GUESS_COUNTER}')

if __name__ == "__main__":
    main()
    