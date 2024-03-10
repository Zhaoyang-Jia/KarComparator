import numpy as np


def sign(x):
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def reverse_dict(input_dict):
    """
    requires mapping to be bijective
    :param input_dict:
    :return:
    """
    output_dict = {}
    for position, name in input_dict.items():
        output_dict[name] = position
    return output_dict


def generate_parabola(start_x, end_x, peak_x, peak_y, num_points=100):
    # Solve for the coefficients of the parabola
    A = np.array([[start_x ** 2, start_x, 1],
                  [peak_x ** 2, peak_x, 1],
                  [end_x ** 2, end_x, 1]])
    b = np.array([0.5, peak_y, 0.5])
    coef = np.linalg.solve(A, b)

    # Generate x values
    x = np.linspace(start_x, end_x, num_points)
    # Calculate y values for the parabola
    y = coef[0] * x ** 2 + coef[1] * x + coef[2]
    # Return x, y as tuples with precision to the sixth decimal point
    return np.array([(round(x_val, 6), round(y_val, 6)) for x_val, y_val in zip(x, y)])


def generate_circle(start_x, diameter, num_points=100):
    # Calculate the radius from the diameter
    radius = diameter / 2

    # Calculate the center of the circle
    center_x = start_x + radius
    center_y = 0.5  # Fixed y-coordinate

    # Generate angles for the circle
    angles = np.linspace(0, 2 * np.pi, num_points)

    # Generate x and y coordinates for the circle
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)
    # Combine x and y coordinates into points
    points = np.array([[round(x_val, 6), round(y_val, 6)] for x_val, y_val in zip(x, y)])

    # Make sure the first and last points are equal to the starting point
    if points[0][0] != start_x or points[-1][0] != start_x:
        raise RuntimeError('circle coordinate formation error')

    return points
