import numpy as np

# Constants from the user's file
ORBIT_DURATION = 35.0
ANIMATION_RUN_TIME = 5.0
SOLAR_SYSTEM_RADIUS = 3.0
DEG_TO_RAD = np.pi / 180.0
RAD_TO_DEG = 180.0 / np.pi


def get_rotation_matrix_x(theta_deg):
    theta = theta_deg * DEG_TO_RAD
    # Manim's rotation matrix for X axis
    # [1, 0, 0]
    # [0, cos, -sin]
    # [0, sin, cos]
    # Note: User code uses this matrix.
    # When applied to a vector v, output is M . v
    return np.array(
        [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]
    )


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def get_alignment_error(tilt_deg, start_orbit_deg, spin_deg):
    # 1. Setup Constants and Deltas
    orbit_delta_deg = (ANIMATION_RUN_TIME / ORBIT_DURATION) * 360.0
    end_orbit_deg = start_orbit_deg + orbit_delta_deg

    # 2. Define Geometries
    # Our coordinate system:
    # Sun at (0,0,0)
    # Earth Orbit: Originally in XY plane, then Rotated by Tilt around X.
    # Earth Axis: Always (0,0,1) (Global Z) because of how user set up coordinates.
    # (The user rotates the Orbit Plane, effectively keeping Earth upright relative to Camera/Global Z).

    rot_matrix = get_rotation_matrix_x(tilt_deg)

    def get_earth_pos(angle_deg):
        angle = angle_deg * DEG_TO_RAD
        # Flat orbit position
        flat_pos = np.array(
            [
                SOLAR_SYSTEM_RADIUS * np.cos(angle),
                SOLAR_SYSTEM_RADIUS * np.sin(angle),
                0.0,
            ]
        )
        # Tilted position
        # M . v
        return np.dot(rot_matrix, flat_pos)

    # 3. Start State
    # Earth Position
    P_start = get_earth_pos(start_orbit_deg)

    # Sun Vector (Earth -> Sun)
    V_sun_start = -P_start  # (0,0,0) - P_start

    # Observer Vector (attached to Earth)
    # At Start, we assume Observer is facing the Sun.
    # Since Earth Axis is Z, we project V_sun into XY plane to match Longitude.
    # (Solar Day is about Meridian crossing).
    V_sun_start_xy = np.array([V_sun_start[0], V_sun_start[1], 0.0])
    Obs_start = normalize(V_sun_start_xy)

    # 4. End State
    # Move Earth
    P_end = get_earth_pos(end_orbit_deg)

    # Sun Vector
    V_sun_end = -P_end
    V_sun_end_xy = np.array([V_sun_end[0], V_sun_end[1], 0.0])
    Target_Vector = normalize(V_sun_end_xy)

    # Move Observer (Spin)
    # Observer rotates around Z axis (Earth Axis) by spin_deg
    spin_rad = spin_deg * DEG_TO_RAD
    cos_spin = np.cos(spin_rad)
    sin_spin = np.sin(spin_rad)

    # Z-rotation matrix
    spin_matrix = np.array(
        [[cos_spin, -sin_spin, 0], [sin_spin, cos_spin, 0], [0, 0, 1]]
    )

    Obs_end = np.dot(spin_matrix, Obs_start)

    # 5. Measure Error
    # Angle between Obs_end and Target_Vector (both in XY plane, normalized)
    # Dot product
    dot = np.dot(Obs_end, Target_Vector)
    # Clamp for numerical stability
    dot = max(min(dot, 1.0), -1.0)

    angle_error_rad = np.arccos(dot)
    angle_error_deg = angle_error_rad * RAD_TO_DEG

    # Determine sign of error (did we overshoot or undershoot?)
    # Cross product Z component
    cross_z = np.cross(Obs_end, Target_Vector)[2]
    # If cross_z > 0, Target is "ahead" (counter-clockwise) of Observer -> We reflected/undershot?
    # Wait, usually spin is CCW. If Target is CCW of Observer, we need MORE spin -> Undershot.

    direction = (
        "Undershot (Need more spin)" if cross_z > 0 else "Overshot (Too much spin)"
    )

    return angle_error_deg, direction, orbit_delta_deg


def get_error_for_duration(duration, tilt, start_orbit_deg):
    # Calculate angles based on duration
    # Spin Rate: 360 deg / ROTATION_DURATION
    # Orbit Rate: 360 deg / ORBIT_DURATION

    spin_rate = 360.0 / 5.0  # ROTATION_DURATION
    orbit_rate = 360.0 / 35.0  # ORBIT_DURATION

    current_spin = spin_rate * duration
    # Orbit delta is how much we moved along the orbit in this time
    orbit_delta = orbit_rate * duration

    err, _, _ = get_alignment_error_raw(
        tilt, start_orbit_deg, current_spin, orbit_delta
    )
    return err, current_spin


def get_alignment_error_raw(tilt_deg, start_orbit_deg, spin_deg, orbit_delta_deg):
    # Helper to avoid recalculating rates in the optimization loop
    # 2. Define Geometries
    end_orbit_deg = start_orbit_deg + orbit_delta_deg
    rot_matrix = get_rotation_matrix_x(tilt_deg)

    def get_earth_pos(angle_deg):
        angle = angle_deg * DEG_TO_RAD
        flat_pos = np.array(
            [
                SOLAR_SYSTEM_RADIUS * np.cos(angle),
                SOLAR_SYSTEM_RADIUS * np.sin(angle),
                0.0,
            ]
        )
        return np.dot(rot_matrix, flat_pos)

    # 3. Start State
    P_start = get_earth_pos(start_orbit_deg)
    V_sun_start = -P_start
    V_sun_start_xy = np.array([V_sun_start[0], V_sun_start[1], 0.0])
    Obs_start = normalize(V_sun_start_xy)

    # 4. End State
    P_end = get_earth_pos(end_orbit_deg)
    V_sun_end = -P_end
    V_sun_end_xy = np.array([V_sun_end[0], V_sun_end[1], 0.0])
    Target_Vector = normalize(V_sun_end_xy)

    # Move Observer (Spin)
    spin_rad = spin_deg * DEG_TO_RAD
    cos_spin = np.cos(spin_rad)
    sin_spin = np.sin(spin_rad)
    spin_matrix = np.array(
        [[cos_spin, -sin_spin, 0], [sin_spin, cos_spin, 0], [0, 0, 1]]
    )
    Obs_end = np.dot(spin_matrix, Obs_start)

    dot = np.dot(Obs_end, Target_Vector)
    dot = max(min(dot, 1.0), -1.0)
    angle_error_rad = np.arccos(dot)
    return (
        angle_error_rad * RAD_TO_DEG,
        Obs_end,
        Target_Vector,
    )  # Reusing return signature slightly differently


# Re-implementing the original get_alignment_error to use the raw helper for backward compatibility/logging
def get_alignment_error(tilt_deg, start_orbit_deg, spin_deg):
    # This was assuming fixed time of 5.0s
    orbit_delta_deg = (5.0 / 35.0) * 360.0
    err, obs_end, target = get_alignment_error_raw(
        tilt_deg, start_orbit_deg, spin_deg, orbit_delta_deg
    )

    # Direction logic
    cross_z = np.cross(obs_end, target)[2]
    direction = "Undershot" if cross_z > 0 else "Overshot"
    return err, direction, orbit_delta_deg


def find_optimal_duration(tilt, start_angle):
    # Theoretical Mean Solar Day Duration
    # w_rot = 2pi/5, w_orb = 2pi/35
    # T = 2pi / (w_rot - w_orb) = 1 / (1/5 - 1/35) = 1 / (6/35) = 35/6 = 5.8333...
    mean_duration = 35.0 / 6.0

    best_t = 0
    min_err = 999.0

    # Search around mean duration +/- 2 seconds
    search_range = np.arange(mean_duration - 1.5, mean_duration + 1.5, 0.01)

    for t in search_range:
        err, _ = get_error_for_duration(t, tilt, start_angle)
        if err < min_err:
            min_err = err
            best_t = t

    # Fine tune
    for t in np.arange(best_t - 0.02, best_t + 0.02, 0.0001):
        err, _ = get_error_for_duration(t, tilt, start_angle)
        if err < min_err:
            min_err = err
            best_t = t

    return best_t, min_err


def run_tests():
    # Mean Solar Day Duration: ~5.833s
    mean_duration = 35.0 / 6.0

    tests = [
        {"tilt": 70, "start": 150},
        {"tilt": 70, "start": 60},
        {"tilt": 23.5, "start": 150},
        {"tilt": 23.5, "start": 60},
        {"tilt": 40.5, "start": 60},
        {"tilt": 0, "start": 0},
    ]

    print(
        f"{'Tilt':<5} | {'Start':<5} | {'Duration (t)':<12} | {'Spin (deg)':<12} | {'Diff from Mean (s)':<20}"
    )
    print("-" * 80)

    for t in tests:
        opt_duration, opt_err = find_optimal_duration(t["tilt"], t["start"])

        # Calculate resulting spin
        spin_rate = 360.0 / 5.0
        final_spin = spin_rate * opt_duration

        diff = opt_duration - mean_duration

        print(
            f"{t['tilt']:<5} | {t['start']:<5} | {opt_duration:<12.5f} | {final_spin:<12.5f} | {diff:+.5f}"
        )


if __name__ == "__main__":
    run_tests()
