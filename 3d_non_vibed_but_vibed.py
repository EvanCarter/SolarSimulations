# Run with: ./.venv/bin/manim -pql 3d_non_vibed_but_vibed.py OrbitRotationTransformation
from manim import *


class SiderealVsSolarNoTilt(ThreeDScene):
    def setup(self):
        # Constants
        self.EARTH_RADIUS = 0.3
        self.ROTATION_DURATION = 5
        self.NUM_SPINS = 1.0
        self.SOLAR_SYSTEM_RADIUS = 3.0
        self.ORBIT_DURATION = 120.0  # Arbitrary slow orbit for now
        self.POINTER_LENGTH = 0.7

        # Camera Setup for 3D
        self.set_camera_orientation(phi=65 * DEGREES, theta=-90 * DEGREES)

        # self.sun = Circle(radius=1).set_color(YELLOW).set_opacity(0.3)
        # turn sun into sphere
        self.sun = Sphere(radius=0.5).set_color(YELLOW).set_opacity(0.3)

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

        # 3. Solar Day
        # Earth has to rotate extra to point back at the sun
        # We calculate the time t where: w_rot * t - w_orb * t = 2 * PI
        w_rot = (2 * PI) / self.ROTATION_DURATION  # angular velocity of spin
        w_orb = (2 * PI) / self.ORBIT_DURATION  # angular velocity of orbit

        # t = 2PI / (w_rot - w_orb)
        solar_day_duration = (2 * PI) / (w_rot - w_orb)

        # Total rotation angle needed = w_rot * t
        solar_day_angle = w_rot * solar_day_duration

        run_day_cycle(
            target_angle=solar_day_angle,
            run_time=solar_day_duration,
            show_solar_arrow=True,
        )


#  now do the scene with the arrow pointing the the sun
class OrbitRotationTransformation(ThreeDScene):
    def construct(self):
        # Define Rotation Angle (23.5 degrees)
        theta = 23.5 * DEGREES

        EARTH_RADIUS = 0.3

        # 1. Setup Camera
        self.set_camera_orientation(phi=75 * DEGREES, theta=10 * DEGREES)

        # Rotation Visualization
        rotation_tracker = ValueTracker(0)

        rotation_label = Tex("Rotation Amount:", font_size=24).to_corner(UR, buff=1)

        number = DecimalNumber(0, font_size=24).next_to(rotation_label, RIGHT)
        self.add_fixed_in_frame_mobjects(rotation_label, number)

        def update_number(m):
            m.set_value(rotation_tracker.get_value() * 180 / PI)
            # when commented this moves it to z=0, +x, +y (aka in the 3d environment)
            # rather than just the top right of the frame
            # seems like a manim bug
            self.add_fixed_in_frame_mobjects(m)

        number.add_updater(update_number)

        # 2. Setup Axes
        axes = ThreeDAxes()
        labels = axes.get_axis_labels(
            Text("x").scale(0.7), Text("y").scale(0.7), Text("z").scale(0.7)
        )
        self.add(axes, labels)

        # 3. Create Orbit on X-Z Plane
        # Initial orbit: flat on the ground (X-Z plane)
        # We start with a circle in X-Y plane and rotate it to X-Z plane
        radius = 3.0
        orbit = Circle(radius=radius, color=BLUE)

        # create sun
        sun = Sphere(radius=0.4).set_color(YELLOW).set_opacity(0.3)

        # 4. Create Earth Group (Sphere + Equator)
        earth_sphere = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E).set_opacity(1)
        # Equator: Circle of same radius, red
        earth_equator = Circle(radius=EARTH_RADIUS, color=RED).set_stroke(width=2)

        # Add a reference line ("stick") to visualize rotation
        # Starts at center, points out to radius + bit more
        stick = Line(
            start=np.array([0, -EARTH_RADIUS, 0]),
            end=np.array([0, -EARTH_RADIUS - 0.7, 0]),
            color=RED,
        ).set_stroke(width=4)

        self.earth_group = Group(earth_sphere, earth_equator, stick)

        # Pre-rotate the Earth group by -theta so that it effectively "starts" tilted relative to the orbit's future frame
        # We rotate about ORIGIN to keep the sphere at (0,0,0) so we can shift it cleanly later
        self.earth_group.rotate(-theta, axis=RIGHT, about_point=ORIGIN)

        # Position at (0, radius, 0)
        # Use shift instead of move_to because move_to centers the bounding box (which is offset due to the stick)
        self.earth_group.shift(np.array([0, radius, 0]))

        # Define Solar Arrow
        self.solar_arrow = Arrow(start=ORIGIN, end=UP, buff=0, color=YELLOW)

        def update_solar_arrow(mob):
            earth_center = earth_sphere.get_center()
            sun_center = sun.get_center()
            vec = sun_center - earth_center

            # FLATTEN THE Z PART SO IT JUST POINTS INWARD RATHER THAN
            # ANGLED APPROPARITATLY TOWARD SUN
            vec[2] = 0
            if np.linalg.norm(vec) > 0.001:
                unit_vec = normalize(vec)
                start = earth_center + unit_vec * EARTH_RADIUS
                end = start + unit_vec * 0.7
                mob.put_start_and_end_on(start, end)

        self.solar_arrow.add_updater(update_solar_arrow)

        # 5. Define Rotation Matrix
        # Rotation around X-axis (fixing typo in comment: matrix is for X-rot, comment said Y/X mixed)
        # The matrix provided:
        # [ 1  0  0 ]
        # [ 0  cos  -sin ]
        # [ 0  sin   cos ]
        # This IS a rotation about the X-axis.

        # Manim's matrix notation is [[row1], [row2], [row3]]
        rotation_matrix = [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]

        self.play(FadeIn(sun), Create(orbit))
        self.play(FadeIn(self.earth_group), FadeIn(self.solar_arrow))
        self.wait(1)

        # 6. Apply Transformation
        # ApplyMatrix applies the linear transformation to the points of the Mobjects
        self.play(
            ApplyMatrix(rotation_matrix, orbit),
            ApplyMatrix(rotation_matrix, self.earth_group),
            run_time=3,
        )

        self.wait(1)

        # 7. Animate Orbit + Rotation
        # We need an orbit tracker to drive the position along the ring
        orbit_tracker = ValueTracker(
            PI / 2
        )  # Start at top (matching the initial placement)

        # Convert our list-of-lists matrix to a numpy array for the updater
        rot_mat_np = np.array(rotation_matrix)

        # Initialize tracking for rotation delta
        self.earth_group.old_angle = rotation_tracker.get_value()

        def move_earth_on_orbit(mob):
            # 1. Calculate target position on the original (flat) circle
            # corresponding to the current orbit_tracker value
            angle = orbit_tracker.get_value()
            flat_pos = np.array([radius * np.cos(angle), radius * np.sin(angle), 0])

            # 2. Apply the same tilt transformation to this position
            # (matrix multiplication)
            tilted_pos = np.dot(rot_mat_np, flat_pos)

            # 3. Move the group to this new position
            # We want the earth_sphere's center to be at 'tilted_pos'.
            # Since 'mob' is the Group (Earth+Stick), we shift the whole group
            # by the difference between where the sphere IS and where it SHOULD BE.
            current_center = earth_sphere.get_center()
            shift_vec = tilted_pos - current_center
            mob.shift(shift_vec)

        def spin_earth(mob):
            # Calculate how much to rotate in this frame
            current_val = rotation_tracker.get_value()
            old_val = mob.old_angle
            delta = current_val - old_val
            mob.old_angle = current_val

            # Rotate around the earth sphere's center (axis=OUT is global Z)
            # Since the sphere moves, we must get its center every frame
            mob.rotate(delta, axis=OUT, about_point=earth_sphere.get_center())

        # Add updaters
        self.earth_group.add_updater(move_earth_on_orbit)
        self.earth_group.add_updater(spin_earth)

        self.play(
            orbit_tracker.animate.increment_value(
                2 * PI / 24
            ),  # 1/24th of an orbit (approx 15 degrees)
            rotation_tracker.animate.increment_value(2 * PI),  # One full rotation
            run_time=5,
            rate_func=linear,
        )

        self.earth_group.remove_updater(move_earth_on_orbit)
        self.earth_group.remove_updater(spin_earth)

        self.wait(2)

        self.move_camera(phi=0, theta=0, run_time=4)
        wait(2)
        self.move_camera(phi=75 * DEGREES, theta=10 * DEGREES, run_time=4)
