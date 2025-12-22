# Run with: ./.venv/bin/manim -pql 3d_non_vibed_but_vibed.py OrbitRotationTransformation
from manim import *


# Constants
EARTH_RADIUS = 0.3
SUN_RADIUS = 0.5
SOLAR_SYSTEM_RADIUS = 3.0
ORBIT_DURATION = 35.0
ROTATION_DURATION = 5.0
POINTER_LENGTH = 0.7

# SIDEREAL LOGIC REMOVED


# 420 for the stationary one

#  WITH A TILT OF 70
# AT A START AT 150 a 378.82 DEGREES GIVES US ONE SOLAR DAY
#  AT A START AT 60 a 468.29 DEGREES GIVES US ONE SOLAR DAY

AXIAL_TILT = 70 * DEGREES
# START_ORBIT_ANGLE = 150 * DEGREES
# TARGET_SPIN_RADIANS = 378.82 * DEGREES
from verify_solar_day import find_optimal_duration

AXIAL_TILT_DEG = 70
START_ORBIT_ANGLE_DEG = 150

AXIAL_TILT = AXIAL_TILT_DEG * DEGREES
START_ORBIT_ANGLE = START_ORBIT_ANGLE_DEG * DEGREES


class SiderealVsSolarNoTilt(ThreeDScene):
    def setup(self):
        # Constants
        self.EARTH_RADIUS = EARTH_RADIUS
        self.ROTATION_DURATION = ROTATION_DURATION
        self.NUM_SPINS = 1.0
        self.SOLAR_SYSTEM_RADIUS = SOLAR_SYSTEM_RADIUS
        self.ORBIT_DURATION = ORBIT_DURATION
        self.POINTER_LENGTH = POINTER_LENGTH

        # Camera Setup for 3D
        self.set_camera_orientation(phi=55 * DEGREES, theta=3 * DEGREES)

        # self.sun = Circle(radius=1).set_color(YELLOW).set_opacity(0.3)
        # turn sun into sphere
        self.sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_opacity(0.3)

        self.planetary_rim = (
            Circle(radius=self.SOLAR_SYSTEM_RADIUS).set_color(GRAY).set_stroke(width=1)
        )
        self.planetary_rim.set_fill(GRAY, opacity=0.2)

        self.earth = Sphere(radius=self.EARTH_RADIUS).set_color(BLUE_E).set_opacity(1)
        earth_center = np.array([0, self.SOLAR_SYSTEM_RADIUS, 0])
        self.earth.move_to(earth_center)

        # Line starts from the edge of the earth (center - radius in y)
        line_start = earth_center - np.array([0, self.EARTH_RADIUS, 0])
        # Original line was ~0.7 units long (2.5 to 1.8)
        line_end = line_start - np.array([0, self.POINTER_LENGTH, 0])
        self.earth_line = Line(line_start, line_end)
        self.earth_line.set_color(YELLOW_A)

        self.earth_and_stick = Group(self.earth, self.earth_line)

        # Initialize orbit angle to PI/2 (starting position (0, 3, 0))
        self.orbit_tracker = ValueTracker(PI / 2)
        self.rotation_tracker = ValueTracker(0)

        # Initialize tracking attribute for rotation delta
        self.earth_and_stick.old_rotation_angle = 0

        self.arrow_pointing_to_sun = Arrow(start=ORIGIN, end=UP, buff=0).set_color(
            YELLOW_A
        )

        # Spin Label
        self.spin_label = (
            Tex("Total Amount Spun: ", font_size=24).to_corner(UR).shift(LEFT * 1.5)
        )
        self.spin_value = DecimalNumber(0, num_decimal_places=1, font_size=24).next_to(
            self.spin_label, RIGHT
        )

        self.add_fixed_in_frame_mobjects(self.spin_label, self.spin_value)

    def construct(self):
        ###### UPDATERS ######
        def orbit_around_sun(obj: Mobject, dt):
            # Update angle
            orbit_rate = (2 * PI) / self.ORBIT_DURATION
            self.orbit_tracker.increment_value(orbit_rate * dt)
            angle = self.orbit_tracker.get_value()

            # Calculate target position for the Earth's center
            target_pos = self.sun.get_center() + np.array(
                [
                    self.SOLAR_SYSTEM_RADIUS * np.cos(angle),
                    self.SOLAR_SYSTEM_RADIUS * np.sin(angle),
                    0,
                ]
            )

            # Shift the entire group so that 'earth' moves to 'target_pos'
            # We use shift instead of move_to because move_to would center the VGroup,
            # which is not the same as centering the Earth due to the line.
            shift_vec = target_pos - self.earth.get_center()
            obj.shift(shift_vec)

        def rotate_earth(obj: Mobject):
            current_angle = self.rotation_tracker.get_value()
            # Calculate delta from previous frame
            # We use getattr/setattr or the attribute we initialized
            old_angle = getattr(obj, "old_rotation_angle", 0)
            delta = current_angle - old_angle

            obj.rotate(
                delta,
                about_point=self.earth.get_center(),
            )
            # set current to be old
            obj.old_rotation_angle = current_angle

        def update_solar_arrow(mob: Arrow, dt=0):
            earth_ctr = self.earth.get_center()
            sun_ctr = self.sun.get_center()
            # Vector from Earth to Sun
            vec = sun_ctr - earth_ctr
            unit_vec = normalize(vec)

            # Emerge from Earth's crust
            start = earth_ctr + unit_vec * self.EARTH_RADIUS
            # Same length as the sidereal line
            end = start + unit_vec * self.POINTER_LENGTH

            mob.put_start_and_end_on(start, end)

        def update_spin_label(mob):
            # Convert radians to degrees
            val = self.rotation_tracker.get_value() * 180 / PI
            mob.set_value(val)
            self.add_fixed_in_frame_mobjects(mob)

        #### END UPDATERS ####

        # Create
        # self.earth_and_stick.save_state() # Sphere doesn't support save_state/restore well
        # Initialize attribute
        self.earth_and_stick.old_rotation_angle = 0

        self.play(Create(self.sun))
        self.play(
            Create(self.planetary_rim),
            FadeIn(self.earth_and_stick),
        )
        self.add(self.spin_label, self.spin_value)

        def run_day_cycle(
            target_angle: float,
            run_time: float = self.ROTATION_DURATION,
            show_solar_arrow: bool = False,
        ) -> None:
            # Resets
            # 1. Reset Position (Orbit to PI/2)
            start_orbit_angle = PI / 2
            target_pos = self.sun.get_center() + np.array(
                [
                    self.SOLAR_SYSTEM_RADIUS * np.cos(start_orbit_angle),
                    self.SOLAR_SYSTEM_RADIUS * np.sin(start_orbit_angle),
                    0,
                ]
            )
            dist_to_reset = target_pos - self.earth.get_center()
            self.earth_and_stick.shift(dist_to_reset)

            # 2. Reset Rotation
            current_rot_accum = getattr(self.earth_and_stick, "old_rotation_angle", 0)
            self.earth_and_stick.rotate(
                -current_rot_accum,
                about_point=self.earth.get_center(),  # Rotate back around center
            )

            self.orbit_tracker.set_value(start_orbit_angle)
            self.rotation_tracker.set_value(0)
            self.earth_and_stick.old_rotation_angle = 0

            # Handle Arrow
            if show_solar_arrow:
                update_solar_arrow(self.arrow_pointing_to_sun)
                self.arrow_pointing_to_sun.add_updater(update_solar_arrow)
                self.add(self.arrow_pointing_to_sun)
            else:
                self.arrow_pointing_to_sun.remove_updater(update_solar_arrow)
                self.remove(self.arrow_pointing_to_sun)

            # Add updaters
            self.earth_and_stick.add_updater(rotate_earth)
            self.earth_and_stick.add_updater(orbit_around_sun)
            self.spin_value.add_updater(update_spin_label)

            # Run Animation
            self.play(
                self.rotation_tracker.animate.set_value(target_angle),
                run_time=run_time,
                rate_func=linear,
            )

            # Clean up
            self.earth_and_stick.remove_updater(rotate_earth)
            self.earth_and_stick.remove_updater(orbit_around_sun)
            self.spin_value.remove_updater(update_spin_label)
            if show_solar_arrow:
                self.arrow_pointing_to_sun.remove_updater(update_solar_arrow)

            self.wait(2)

        # 1. Sidereal Day (360 degrees, points to distinct stars, not sun)
        run_day_cycle(target_angle=2 * PI * self.NUM_SPINS, show_solar_arrow=False)

        # 2. Sidereal Day with Reference Arrow
        # Shows that 360 rotation doesn't arrive back at pointing to sun
        run_day_cycle(target_angle=2 * PI * self.NUM_SPINS, show_solar_arrow=True)

        run_day_cycle(
            target_angle=SOLAR_DAY_ANGLE_SIDEREAL,
            run_time=SOLAR_DAY_DURATION_SIDEREAL,
            show_solar_arrow=True,
        )


#  now do the scene with the arrow pointing the the sun
class OrbitRotationTransformation(ThreeDScene):

    def camera_to_origin(self):
        self.move_camera(phi=0, theta=0, focal_distance=400, run_time=2.5)
        # self.move_camera(phi=0, theta=0, focal_distance=1000, run_time=1.5)
        self.wait(0.5)
        self.move_camera(
            phi=self.camera_phi_angle, theta=self.camera_theta_angle, run_time=1.8
        )
        self.wait(0.5)

    def get_meridian_aligned_to_sun(self, earth_sphere, earth_axis, sun):
        # 1. Get current orientation vectors
        # Vector pointing from Earth center to North Pole
        curr_earth_axis_vec = earth_axis.get_end() - earth_axis.get_start()
        n_axis = normalize(curr_earth_axis_vec)

        # Vector pointing from Earth center to Sun
        vec_to_sun = sun.get_center() - earth_sphere.get_center()

        # 2. Project vec_to_sun onto the plane perpendicular to the axis
        # proj = v - (v . n) * n
        proj_to_sun = vec_to_sun - np.dot(vec_to_sun, n_axis) * n_axis

        # Normalize to get the "forward" direction on the Earth's equator plane
        n_forward = normalize(proj_to_sun)

        # 3. Calculate "Right" vector (East direction from local perspective)
        # Right hand rule: forward x axis ? No, usually we want a right-handed basis.
        # If Z is axis, X is forward, then Y is Z x X ... wait.
        # Let's just define a basis matrix where:
        # X' (mapped from [1,0,0]) -> n_forward (Points to Sun in equatorial plane)
        # Y' (mapped from [0,1,0]) -> n_axis (Points North)
        # Z' (mapped from [0,0,1]) -> n_right (Points East-ish)
        n_right = np.cross(n_forward, n_axis)

        # Create the transformation matrix
        # M = [col1, col2, col3]
        basis_matrix = np.array([n_forward, n_axis, n_right]).T

        # 4. Create Meridian
        # A semi-circle going from South Pole to North Pole? Or North to South?
        # Standard Arc starts at 0 rad (X-axis).
        # range [-PI/2, PI/2] would be right half of circle (4th and 1st quadrants).
        # We want it to pass through the "forward" direction.
        # If we map Standard X to Forward, Standard Y to Axis.
        # Arc from -PI/2 to PI/2 in standard XY plane goes from (0,-R) to (0,R) through (R,0).
        # This corresponds perfectly to South Pole -> North Pole through Forward point.
        meridian = Arc(
            radius=EARTH_RADIUS,
            start_angle=-PI / 2,
            angle=PI,
            color=YELLOW,
        ).set_stroke(width=4)

        # Apply the basis transformation
        # We can use apply_matrix
        meridian.apply_matrix(basis_matrix)

        # Move to Earth center (apply_matrix rotates around origin, assuming arc was created at origin)
        meridian.shift(earth_sphere.get_center())

        return meridian

    def run_orbit_scenario(
        self,
        *,
        axial_tilt_deg,
        start_orbit_angle_deg,
        is_flat_orbit,
        orbit_mobject,
        sun_mobject,
        show_earth_angle_rotation: bool,
    ):
        """
        Runs a single orbit scenario.
        orbit_mobject is passed in so it can be reused (and tilted/untilted) without being destroyed.
        """
        AXIAL_TILT = axial_tilt_deg * DEGREES
        START_ORBIT_ANGLE = start_orbit_angle_deg * DEGREES

        # Calculate dynamic duration
        try:
            scenario_solar_day_duration, _ = find_optimal_duration(
                axial_tilt_deg, start_orbit_angle_deg
            )
            print(f"Scenario Duration: {scenario_solar_day_duration}")
        except Exception as e:
            print(f"Error calculating duration: {e}")
            scenario_solar_day_duration = 5.0

        # Derived Physics
        w_rot = (2 * PI) / ROTATION_DURATION
        TARGET_SPIN_RADIANS = w_rot * scenario_solar_day_duration
        ANIMATION_RUN_TIME = scenario_solar_day_duration

        # Scene Group for easy cleanup
        scene_group = Group()

        # Rotation Visualization (UI)
        rotation_tracker = ValueTracker(0)
        rotation_label = Tex("Rotation Amount:", font_size=24).to_corner(UR, buff=1)
        number = DecimalNumber(0, font_size=24).next_to(rotation_label, RIGHT)
        # We handle these separately for cleanup since they are fixed in frame

        def update_number(m):
            m.set_value(rotation_tracker.get_value() * 180 / PI)
            self.add_fixed_in_frame_mobjects(m)

        number.add_updater(update_number)

        # Create Earth Group (Sphere + Equator)
        earth_sphere = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E).set_opacity(1)
        earth_equator = Circle(radius=EARTH_RADIUS, color=RED).set_stroke(width=2)

        EARTH_START_X = SOLAR_SYSTEM_RADIUS * np.cos(START_ORBIT_ANGLE)
        EARTH_START_Y = SOLAR_SYSTEM_RADIUS * np.sin(START_ORBIT_ANGLE)

        earth_axis = Line(
            start=np.array([0, 0, -EARTH_RADIUS * 1.4]),
            end=np.array([0, 0, EARTH_RADIUS * 1.4]),
            color=GREEN,
        ).set_stroke(width=2)

        earth_group = Group(earth_sphere, earth_equator, earth_axis)
        scene_group.add(earth_group)

        # Position Earth
        earth_group.shift(np.array([EARTH_START_X, EARTH_START_Y, 0]))

        # Solar Arrow
        solar_arrow = Arrow(start=ORIGIN, end=UP, buff=0, color=YELLOW)
        scene_group.add(solar_arrow)

        def update_solar_arrow(mob):
            earth_center = earth_sphere.get_center()
            sun_center = sun_mobject.get_center()
            vec = sun_center - earth_center

            # Flatten Z for visual clarity
            vec[2] = 0
            if np.linalg.norm(vec) > 0.001:
                unit_vec = normalize(vec)
                start = earth_center + unit_vec * EARTH_RADIUS
                end = start + unit_vec * 0.7
                mob.put_start_and_end_on(start, end)

        solar_arrow.add_updater(update_solar_arrow)

        # Define Rotation Matrix
        theta = AXIAL_TILT
        rotation_matrix = [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]

        self.add_fixed_in_frame_mobjects(rotation_label, number)

        # Animation Start
        if show_earth_angle_rotation:
            self.play(FadeIn(earth_group), FadeIn(solar_arrow))
            self.wait(1)
            # Tilt the earth
            self.play(
                earth_group.animate.rotate(
                    -theta,
                    axis=RIGHT,
                    about_point=earth_sphere.get_center(),
                ),
                run_time=3,
            )
        else:
            earth_group.rotate(
                -theta, axis=RIGHT, about_point=earth_sphere.get_center()
            )

        meridian = self.get_meridian_aligned_to_sun(
            earth_sphere, earth_axis, sun_mobject
        )
        earth_group.add(meridian)
        scene_group.add(meridian)
        if show_earth_angle_rotation:
            self.play(FadeIn(meridian), run_time=2)
        else:
            self.play(FadeIn(earth_group), FadeIn(solar_arrow))

        # Matrix Transformation
        # If not a flat orbit, we tilt the entire coordinate system (Orbit + Earth Group)
        if not is_flat_orbit:
            self.play(
                ApplyMatrix(rotation_matrix, orbit_mobject),
                ApplyMatrix(rotation_matrix, earth_group),
                run_time=3,
            )

        # Camera Move
        self.camera_to_origin()

        # Ambient Rotation
        self.begin_ambient_camera_rotation(rate=-0.01)

        # Orbit + Spin Animation
        orbit_tracker = ValueTracker(START_ORBIT_ANGLE)
        rot_mat_np = np.array(rotation_matrix)
        earth_group.old_angle = rotation_tracker.get_value()

        def move_earth_on_orbit(mob):
            angle = orbit_tracker.get_value()
            flat_pos = np.array(
                [
                    SOLAR_SYSTEM_RADIUS * np.cos(angle),
                    SOLAR_SYSTEM_RADIUS * np.sin(angle),
                    0,
                ]
            )
            current_center = earth_sphere.get_center()

            if is_flat_orbit:
                shift_vec = flat_pos - current_center
            else:
                tilted_pos = np.dot(rot_mat_np, flat_pos)
                shift_vec = tilted_pos - current_center

            mob.shift(shift_vec)

        def spin_earth(mob):
            current_val = rotation_tracker.get_value()
            old_val = mob.old_angle
            delta = current_val - old_val
            mob.old_angle = current_val

            # Axis of rotation
            if is_flat_orbit:
                # Tilt is local to Earth, but orbit is flat.
                # Earth was rotated by -theta around RIGHT.
                # So its axis is [0, sin(theta), cos(theta)]
                axis = np.array([0, np.sin(theta), np.cos(theta)])
            else:
                axis = OUT

            mob.rotate(delta, axis=axis, about_point=earth_sphere.get_center())

        earth_group.add_updater(move_earth_on_orbit)
        earth_group.add_updater(spin_earth)

        # Calculate increment
        orbit_increment = (ANIMATION_RUN_TIME / ORBIT_DURATION) * TAU

        self.play(
            orbit_tracker.animate.increment_value(orbit_increment),
            rotation_tracker.animate.increment_value(TARGET_SPIN_RADIANS),
            run_time=ANIMATION_RUN_TIME,
            rate_func=linear,
        )

        # Cleanup Updaters
        self.stop_ambient_camera_rotation()
        earth_group.remove_updater(move_earth_on_orbit)
        earth_group.remove_updater(spin_earth)
        solar_arrow.remove_updater(update_solar_arrow)
        number.remove_updater(update_number)

        self.camera_to_origin()

        # FADE OUT
        # Remove fixed elements

        # Restore Orbit if needed (Inverse Pivot)
        cleanup_anims = [FadeOut(scene_group), FadeOut(number), FadeOut(rotation_label)]

        self.play(*cleanup_anims, run_time=1.5)
        self.remove_fixed_in_frame_mobjects(rotation_label, number)
        if not is_flat_orbit:
            # We need to un-tilt the orbit
            # Inverse of rotation_matrix (it's orthogonal, so transpose)
            inv_matrix = np.array(rotation_matrix).T
            # Applying matrix to Mobject
            self.play(ApplyMatrix(inv_matrix, orbit_mobject), run_time=1.5)
        self.wait(0.5)

    def construct(self):
        self.camera_phi_angle = 55 * DEGREES
        self.camera_theta_angle = 3 * DEGREES

        # Setup Global Persistent Objects
        self.set_camera_orientation(
            phi=self.camera_phi_angle, theta=self.camera_theta_angle
        )

        axes = ThreeDAxes()
        labels = axes.get_axis_labels(
            Text("x").scale(0.7), Text("y").scale(0.7), Text("z").scale(0.7)
        )
        self.add(axes, labels)

        orbit = Circle(radius=SOLAR_SYSTEM_RADIUS, color=BLUE)
        sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_opacity(0.3)
        self.add(sun)
        self.add(orbit)
        self.play(FadeIn(sun), FadeIn(orbit), run_time=1)
        self.wait(0.5)

        # Scenario 1: Tilted Orbit at 70 degrees
        # self.run_orbit_scenario(
        #     axial_tilt_deg=70,
        #     start_orbit_angle_deg=60,
        #     is_flat_orbit=True,
        #     orbit_mobject=orbit,
        #     sun_mobject=sun,
        #     show_earth_angle_rotation=True,
        # )

        self.run_orbit_scenario(
            axial_tilt_deg=70,
            start_orbit_angle_deg=60,
            is_flat_orbit=False,
            orbit_mobject=orbit,
            sun_mobject=sun,
            show_earth_angle_rotation=False,
        )
