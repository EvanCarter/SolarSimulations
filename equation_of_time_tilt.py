from manim import *
import numpy as np


class EquationOfTimeTilt(ThreeDScene):
    def construct(self):
        # 1. Setup Scene
        # Tweak camera to see both sides better or just rely on shifting objects
        # self.set_camera_orientation(phi=60 * DEGREES, theta=-30 * DEGREES)
        # Zoom out slightly by default? Or just position things well.
        self.set_camera_orientation(phi=60 * DEGREES, theta=-30 * DEGREES)

        # Constants
        TILT = 23.44 * DEGREES
        ORBIT_RADIUS = 2.5  # Slightly smaller to fit in half
        SUN_RADIUS = 0.4
        EARTH_RADIUS = 0.15

        # 2. Left Side: The Orbit System
        orbit_center = LEFT * 3.5

        sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_sheen(0.5, "ul")
        sun.move_to(orbit_center)
        sun_glow = (
            Sphere(radius=SUN_RADIUS * 1.5)
            .set_color(YELLOW)
            .set_opacity(0.2)
            .move_to(orbit_center)
        )

        # Ecliptic (White)
        ecliptic = Circle(radius=ORBIT_RADIUS).set_color(WHITE).set_opacity(0.3)
        ecliptic.move_to(orbit_center)
        ecliptic_line = Circle(radius=ORBIT_RADIUS).set_color(WHITE).set_stroke(width=1)
        ecliptic_line.move_to(orbit_center)

        # Celestial Equator (Blue, Tilted)
        celestial_equator = Circle(radius=ORBIT_RADIUS).set_color(BLUE)
        celestial_equator.rotate(TILT, axis=RIGHT)
        celestial_equator.move_to(orbit_center)

        # Labels for Seasons
        # Vernal (0), Summer (PI/2), Autumn (PI), Winter (3PI/2)
        # Positions are relative to orbit_center
        # Note: Ecliptic is flat in XY.
        # Check angle alignment. At 0, we want Vernal Equinox.
        label_dist = ORBIT_RADIUS + 0.5
        vernal_pos = orbit_center + RIGHT * label_dist
        summer_pos = orbit_center + UP * label_dist
        autumn_pos = orbit_center + LEFT * label_dist
        winter_pos = orbit_center + DOWN * label_dist

        # We use Text (2D) but place them in 3D. They might not face camera perfectly.
        # billboard() or keeping them fixed_in_frame is better but tracking position is hard.
        # Let's rotate them to lie flat or face camera.
        l_vernal = (
            Text("Vernal\nEquinox", font_size=16)
            .move_to(vernal_pos)
            .rotate(PI / 2, axis=RIGHT)
        )
        l_summer = (
            Text("Summer\nSolstice", font_size=16)
            .move_to(summer_pos)
            .rotate(PI / 2, axis=RIGHT)
        )
        l_autumn = (
            Text("Autumnal\nEquinox", font_size=16)
            .move_to(autumn_pos)
            .rotate(PI / 2, axis=RIGHT)
        )
        l_winter = (
            Text("Winter\nSolstice", font_size=16)
            .move_to(winter_pos)
            .rotate(PI / 2, axis=RIGHT)
        )

        orbit_labels = VGroup(l_vernal, l_summer, l_autumn, l_winter)

        # Earth
        earth = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E)

        # "Red Circle signifying tilt" on the Earth
        # This represents Earth's Equator, which is tilted relative to orbit.
        # Earth is at orbit_center + R.
        # Earth Obliquity is 23.5.
        # If we just rotate the Earth sphere, it's hard to see.
        # Add a ring around Earth.
        earth_equator_ring = (
            Circle(radius=EARTH_RADIUS + 0.02).set_color(RED).set_stroke(width=2)
        )
        # Tilt this ring.
        earth_equator_ring.rotate(TILT + PI / 2, axis=RIGHT)
        # Wait, Earth's axis is tilted 23.5 deg from vertical.
        # So Earth's Equator is tilted 23.5 deg from Horizontal.
        # TILT constant is 23.5 deg.
        # So rotate around Right axis by TILT.
        earth_equator_ring.rotate(TILT, axis=RIGHT)

        earth_group = VGroup(earth, earth_equator_ring)

        # Grouping
        # Note: celestial_equator is already moved.
        orbit_system = VGroup(
            ecliptic, ecliptic_line, celestial_equator, sun, sun_glow, orbit_labels
        )
        self.add(orbit_system, earth_group)

        # 3. Dynamic Elements
        lambda_tracker = ValueTracker(0)

        def get_ecliptic_pos_local(lmbda):
            return ORBIT_RADIUS * np.array([np.cos(lmbda), np.sin(lmbda), 0])

        def get_ecliptic_pos_global(lmbda):
            return orbit_center + get_ecliptic_pos_local(lmbda)

        def get_equator_projection_global(lmbda):
            # Same math, just offset by orbit_center
            alpha = np.arctan2(np.cos(TILT) * np.sin(lmbda), np.cos(lmbda))

            raw_pos = ORBIT_RADIUS * np.array([np.cos(alpha), np.sin(alpha), 0])
            rot_matrix = rotation_matrix(TILT, RIGHT)
            final_pos = np.dot(rot_matrix, raw_pos)
            return orbit_center + final_pos, alpha

        # Updating Earth
        earth_group.add_updater(
            lambda m: m.move_to(get_ecliptic_pos_global(lambda_tracker.get_value()))
        )

        earth_group.to_edge(LEFT, buff=1.0)

        # Projection Line (Yellow dashed)
        projection_line = always_redraw(
            lambda: Line(
                orbit_center,
                get_equator_projection_global(lambda_tracker.get_value())[0],
                color=BLUE_A,
                stroke_width=2,
            )
        )

        # Connection Line (The Difference)
        connection_line = always_redraw(
            lambda: Line(
                earth_group.get_center(),
                get_equator_projection_global(lambda_tracker.get_value())[0],
                color=RED,
                stroke_opacity=0.8,
            )
        )

        self.add(projection_line, connection_line)

        # 4. Right Side: The Graph
        # Position axes on the right side of the screen

        axes = Axes(
            x_range=[0, 2 * PI, PI / 2],
            y_range=[0.8, 1.2, 0.1],
            x_length=5,
            y_length=4,
            axis_config={"color": GREY},
        )

        # Add labels manually using Text
        x_label = Text("Orbit Angle", font_size=18).next_to(axes.x_axis, DOWN)
        y_label = (
            Text("Day Scale", font_size=18).rotate(PI / 2).next_to(axes.y_axis, LEFT)
        )

        graph_group = VGroup(axes, x_label, y_label)
        graph_group.to_edge(RIGHT, buff=1.0)
        graph_group.scale(0.9)

        # Background for graph to verify position
        bg_rect = BackgroundRectangle(
            graph_group, color=BLACK, fill_opacity=0.3, buff=0.2
        )

        # Plot function
        def day_len_factor(lmbda):
            A = np.cos(TILT)
            val = A / (np.cos(lmbda) ** 2 + A**2 * np.sin(lmbda) ** 2)
            return val

        graph = axes.plot(day_len_factor, color=GREEN)

        # Moving Dot
        graph_dot = always_redraw(
            lambda: Dot(
                axes.c2p(
                    lambda_tracker.get_value() % (2 * PI),
                    day_len_factor(lambda_tracker.get_value()),
                ),
                color=RED,
            )
        )

        # Vertical Line on graph
        v_line = always_redraw(
            lambda: axes.get_vertical_line(
                axes.c2p(
                    lambda_tracker.get_value() % (2 * PI),
                    day_len_factor(lambda_tracker.get_value()),
                )
            )
        )

        self.add_fixed_in_frame_mobjects(bg_rect, graph_group, graph, graph_dot, v_line)

        # 5. Animation
        self.wait(1)

        # Slow down: 24 seconds for full year
        self.play(
            lambda_tracker.animate.set_value(2 * PI), run_time=24, rate_func=linear
        )
        self.wait()


class SiderealVsSolarNoTilt(Scene):
    def construct(self):
        # Constants
        ORBIT_RADIUS = 3.0
        SUN_RADIUS = 0.4
        EARTH_RADIUS = 0.3
        # 1 day = 1/36 rotation -> 10 degrees orbit.
        # User wants "one day be 1/36th of a rotation" -> Assuming 1/36 of Orbit for a full day?
        # Actually standard is 365 days. 1/36 is exaggerated.
        ORBIT_ARC_PER_DAY = (360 / 36) * DEGREES  # 10 degrees

        # 1. Setup Objects
        sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_sheen(0.5, "ul")
        sun.move_to(ORIGIN)

        orbit_path = Circle(radius=ORBIT_RADIUS).set_color(WHITE).set_opacity(0.3)
        orbit_path_stroke = Circle(radius=ORBIT_RADIUS, color=WHITE, stroke_width=1)

        # Earth Components (Created once, updated via updater)
        earth_sphere = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E)

        # Lines on Earth
        # Meridian (Red)
        meridian_line = Line(
            DOWN * (EARTH_RADIUS + 0.1),
            UP * (EARTH_RADIUS + 0.3),
            color=RED,
            stroke_width=4,
        )
        # Longitudinal lines
        long_lines = VGroup()
        for i in range(8):
            angle = i * (2 * PI / 8)
            l = Line(
                DOWN * EARTH_RADIUS, UP * EARTH_RADIUS, color=BLUE_A, stroke_width=1
            )
            l.rotate(
                angle, axis=OUT, about_point=ORIGIN
            )  # Rotate around Z (2D plane perpendicular)
            long_lines.add(l)

        # Combine into a clean group centered at local origin
        # Note: We keep them at ORIGIN mostly, then move them.
        earth_clean_group = VGroup(earth_sphere, long_lines, meridian_line)

        # Trackers
        # Orbit Angle: Starts at 0 (Right).
        theta_tracker = ValueTracker(0)
        # Spin Angle: Starts at 0.
        spin_tracker = ValueTracker(0)

        # Updater Function
        def update_earth_position(m):
            theta = theta_tracker.get_value()
            spin = spin_tracker.get_value()

            # 1. Reset to base state (Center at origin, no rotation)
            m.restore()

            # 2. Apply Spin
            # Rotate around ITSELF (Z axis)
            m.rotate(spin, axis=OUT)

            # 3. Move to Orbit Position
            # Position: R * [cos(theta), sin(theta), 0]
            pos = ORBIT_RADIUS * np.array([np.cos(theta), np.sin(theta), 0])
            m.move_to(pos)

            # 4. Rotate Meridian Correction (Initial orientation)
            # We want Meridian to point to Sun (Left) at Theta=0 (Right).
            # Currently Meridian is UP.
            # Rotate 90 deg (PI/2) CCW.
            m.rotate(PI / 2, axis=OUT)

        # Prepare group for restoration
        earth_clean_group.save_state()
        earth_clean_group.add_updater(update_earth_position)

        # Reference Arrows
        # We compute these based on trackers to avoid lag
        def get_earth_center():
            theta = theta_tracker.get_value()
            return ORBIT_RADIUS * np.array([np.cos(theta), np.sin(theta), 0])

        solar_arrow = always_redraw(
            lambda: Arrow(
                start=get_earth_center(),
                end=get_earth_center() + normalize(ORIGIN - get_earth_center()) * 1.0,
                color=YELLOW,
                buff=0,
            )
        )

        sidereal_arrow = always_redraw(
            lambda: Arrow(
                start=get_earth_center(),
                end=get_earth_center() + LEFT * 1.0,  # Fixed star direction: LEFT
                color=PURPLE,
                buff=0,
            )
        )

        # Add to Scene
        self.add(orbit_path, orbit_path_stroke, sun)
        self.add(earth_clean_group)
        self.add(solar_arrow, sidereal_arrow)

        self.wait(2)  # Initial pause

        # Animation Phase 1: Sidereal Day
        # Orbit: 10 degrees
        # Spin: 360 degrees

        self.play(
            theta_tracker.animate.set_value(ORBIT_ARC_PER_DAY),
            spin_tracker.animate.set_value(2 * PI),
            run_time=4,
            rate_func=linear,
        )

        # Pause to show mismatch
        self.next_section("Mismatch")
        self.wait(1)

        # Correct for Solar Day
        # Spin needs to increase by the orbit angle (ORBIT_ARC_PER_DAY)
        # Note: Direction depends on rotation correctness.
        # CCW Orbit + CCW Spin -> Needs EXTRA spin?
        # Let's see: Start at Right. Moved UP-LEFT (CCW).
        # Sun is "Back".
        # Meridian did 360, points Left (Initial Space Dir).
        # Sun is MORE Left (down-left).
        # Need to rotate MORE CCW.

        self.play(spin_tracker.animate.increment_value(ORBIT_ARC_PER_DAY), run_time=1.5)

        self.wait(2)


class SiderealVsSolarSolstice(ThreeDScene):
    def construct(self):
        # Constants
        TILT = 45 * DEGREES
        ORBIT_RADIUS = 3.0
        SUN_RADIUS = 0.4
        EARTH_RADIUS = 0.3
        ORBIT_ARC_PER_DAY = (360 / 36) * DEGREES

        # Camera
        self.set_camera_orientation(phi=60 * DEGREES, theta=0 * DEGREES)  # Front view

        # 1. Setup Objects
        sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_sheen(0.5, "ul")
        sun.move_to(ORIGIN)
        ecliptic = Circle(radius=ORBIT_RADIUS, color=WHITE, stroke_width=1).set_opacity(
            0.3
        )
        self.add(ecliptic, sun)

        # Earth Construction
        # We build it at origin, then transform in updater
        sphere = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E)
        ring = Circle(radius=EARTH_RADIUS + 0.05, color=GREEN).set_stroke(width=2)
        # Pointer: From center to FRONT (Z) in local space.
        pointer = Line(ORIGIN, OUT * (EARTH_RADIUS + 0.3), color=RED, stroke_width=4)

        # Helper: Local Group
        local_earth = VGroup(sphere, ring, pointer)

        # Initial Orientation in Local Space:
        # Axis is UP (Y). Equator is in XZ plane.
        # Ring is in XY plane by default. Rotate to XZ.
        ring.rotate(PI / 2, axis=RIGHT)
        # Pointer is along Z.

        # Trackers
        # Solstice: Earth Axis is maximally tilted towards or away from the Sun.
        # If Axis is Y, tilted around X (RIGHT) towards Z.
        # So Axis is in YZ plane.
        # Sun-Earth line is X axis. (Earth at Left/Right).
        # Then Axis (YZ) is perpendicular to Sun-Earth (X). This is Equinox.
        # Solstice is when Sun-Earth line is in the YZ plane.
        # So Earth is at UP (+Y) or DOWN (-Y).
        # Let's start at UP (PI/2).

        theta_tracker = ValueTracker(PI / 2)
        spin_tracker = ValueTracker(0)

        # Re-implement clean updater
        def update_wrapped(m):
            theta = theta_tracker.get_value()
            spin = spin_tracker.get_value()
            m.restore()  # Resets to local_earth's saved state (Axis Y, Equator XZ, Pointer Z)

            # Pre-Spin Orient: Point Z-pointer to Sun direction.
            # Earth starts at UP (+Y). Sun at DOWN (-Y).
            # Pointer is Z.
            # Rotate -90 X to make Z -> -Y.
            m.rotate(-PI / 2, axis=RIGHT)

            # Spin around Y (Axis)
            m.rotate(spin, axis=UP)

            # Tilt (Axis Y tips towards Z)
            # Rotate around X (RIGHT) by TILT.
            m.rotate(TILT, axis=RIGHT)

            # Move
            pos = ORBIT_RADIUS * np.array([np.cos(theta), np.sin(theta), 0])
            m.move_to(pos)

        local_earth.save_state()
        local_earth.add_updater(update_wrapped)
        self.add(local_earth)

        # Solar Arrow
        solar_arrow = always_redraw(
            lambda: Arrow(
                start=local_earth.get_center(),
                end=local_earth.get_center()
                + normalize(ORIGIN - local_earth.get_center()) * 1.5,
                color=YELLOW,
                buff=0,
            )
        )
        self.add(solar_arrow)

        self.wait(1)

        # Anim: Sidereal Day (360 spin) + Orbit
        # Motion: Left (CCW from UP is LEFT).
        # Theta PI/2 -> PI/2 + ORBIT_ARC_PER_DAY.
        self.play(
            theta_tracker.animate.increment_value(ORBIT_ARC_PER_DAY),
            spin_tracker.animate.set_value(2 * PI),
            run_time=4,
        )

        self.next_section("Mismatch Solstice")
        self.wait(1)

        # Visual Check:
        # At Solstice, Solar Day is LONGER.
        # Means we need MORE spin to face Sun.
        # The extra rotation needed is larger than ORBIT_ARC_PER_DAY.

        self.play(
            spin_tracker.animate.increment_value(
                ORBIT_ARC_PER_DAY * 1.15
            ),  # 15% longer
            run_time=2,
        )
        self.wait(2)


class SiderealVsSolarEquinox(ThreeDScene):
    def construct(self):
        # Constants
        TILT = 45 * DEGREES
        ORBIT_RADIUS = 3.0
        SUN_RADIUS = 0.4
        EARTH_RADIUS = 0.3
        ORBIT_ARC_PER_DAY = (360 / 36) * DEGREES

        # Camera: Top/Side hybrid
        self.set_camera_orientation(phi=60 * DEGREES, theta=-90 * DEGREES)

        # 1. Setup
        sun = Sphere(radius=SUN_RADIUS).set_color(YELLOW).set_sheen(0.5, "ul")
        sun.move_to(ORIGIN)
        ecliptic = Circle(radius=ORBIT_RADIUS, color=WHITE, stroke_width=1).set_opacity(
            0.3
        )
        self.add(ecliptic, sun)

        # Earth
        sphere = Sphere(radius=EARTH_RADIUS).set_color(BLUE_E)
        ring = Circle(radius=EARTH_RADIUS + 0.05, color=GREEN).set_stroke(width=2)
        pointer = Line(ORIGIN, OUT * (EARTH_RADIUS + 0.3), color=RED, stroke_width=4)

        local_earth = VGroup(sphere, ring, pointer)

        # Initial Orientation in Local Space:
        # Axis is UP (Y). Equator is in XZ plane.
        # Ring is in XY plane by default. Rotate to XZ.
        ring.rotate(PI / 2, axis=RIGHT)
        # Pointer is along Z.

        # Trackers
        # Equinox: Earth Axis is Perpendicular to Sun Vector.
        # If Axis is Y, tilted around X (RIGHT) towards Z.
        # So Axis is in YZ plane.
        # Sun-Earth line is X axis. (Earth at Left/Right).
        # Then Axis (YZ) is perpendicular to Sun-Earth (X). This is Equinox.
        # Position: Right (0 radians).

        theta_tracker = ValueTracker(0)
        spin_tracker = ValueTracker(0)

        def update_wrapped(m):
            theta = theta_tracker.get_value()
            spin = spin_tracker.get_value()
            m.restore()  # Resets to local_earth's saved state (Axis Y, Equator XZ, Pointer Z)

            # Pre-Spin Orient: Pointer Z. Earth Right. Sun Left (-X).
            # Want Pointer Left (-X).
            # Rotate 90 Y (Z->X). Then 180 around Y to make X -> -X.
            # Total 90 Y.
            m.rotate(PI / 2, axis=UP)

            # Spin around Y (Axis)
            m.rotate(spin, axis=UP)

            # Tilt (Axis Y tips towards Z)
            # Rotate around X (RIGHT) by TILT.
            m.rotate(TILT, axis=RIGHT)

            # Move
            pos = ORBIT_RADIUS * np.array([np.cos(theta), np.sin(theta), 0])
            m.move_to(pos)

        local_earth.save_state()
        local_earth.add_updater(update_wrapped)
        self.add(local_earth)

        solar_arrow = always_redraw(
            lambda: Arrow(
                start=local_earth.get_center(),
                end=local_earth.get_center()
                + normalize(ORIGIN - local_earth.get_center()) * 1.5,
                color=YELLOW,
                buff=0,
            )
        )
        self.add(solar_arrow)

        self.wait(1)

        # Anim: Orbit Upwards (0 -> ORBIT_ARC_PER_DAY)
        self.play(
            theta_tracker.animate.increment_value(ORBIT_ARC_PER_DAY),
            spin_tracker.animate.set_value(2 * PI),
            run_time=4,
        )

        self.next_section("Mismatch Equinox")
        self.wait(1)

        # Equinox Effect: Solar Day is SHORTER (closer to sidereal).
        # The extra rotation needed is smaller than ORBIT_ARC_PER_DAY.

        self.play(
            spin_tracker.animate.increment_value(
                ORBIT_ARC_PER_DAY * 0.85
            ),  # Slower catch up
            run_time=2,
        )
        self.wait(2)
