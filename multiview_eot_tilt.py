from manim import *
import numpy as np


class MultiViewEquationOfTimeTilt(ThreeDScene):
    def setup_constants(self):
        """Initialize all constants for the scene"""
        self.TILT = 23.44 * DEGREES
        self.ORBIT_RADIUS = 2.8
        self.SUN_RADIUS = 0.4
        self.EARTH_RADIUS = 0.18

        # Colors
        self.SUN_COLOR = YELLOW
        self.EARTH_COLOR = BLUE_E
        self.ECLIPTIC_COLOR = WHITE
        self.EQUATOR_COLOR = RED
        self.AXIS_COLOR = YELLOW

        # Camera
        self.PHI = 75 * DEGREES
        self.THETA = -45 * DEGREES

    def create_orbital_system(self, orbit_center):
        """Create sun and orbit path"""
        sun = Sphere(radius=self.SUN_RADIUS).set_color(self.SUN_COLOR)
        sun.move_to(orbit_center)

        sun_glow = Sphere(radius=self.SUN_RADIUS * 1.3)
        sun_glow.set_color(self.SUN_COLOR).set_opacity(0.15)
        sun_glow.move_to(orbit_center)

        ecliptic = Circle(radius=self.ORBIT_RADIUS, color=self.ECLIPTIC_COLOR)
        ecliptic.set_stroke(width=2, opacity=0.5)
        ecliptic.move_to(orbit_center)

        return VGroup(sun, sun_glow, ecliptic)

    def create_earth(self, with_tilt=False):
        """Create Earth with optional tilt indicators"""
        earth_sphere = Sphere(radius=self.EARTH_RADIUS)
        earth_sphere.set_color(self.EARTH_COLOR)

        # Equatorial ring
        equator_ring = Circle(radius=self.EARTH_RADIUS * 1.1)
        equator_ring.set_color(self.EQUATOR_COLOR).set_stroke(width=3)
        equator_ring.rotate(PI / 2, axis=RIGHT)
        if with_tilt:
            equator_ring.rotate(self.TILT, axis=RIGHT)

        # Rotation axis line
        axis_line = Line(
            DOWN * (self.EARTH_RADIUS * 1.3),
            UP * (self.EARTH_RADIUS * 1.3),
            color=self.AXIS_COLOR,
            stroke_width=4,
        )
        if with_tilt:
            axis_line.rotate(self.TILT, axis=RIGHT)

        earth_group = VGroup(earth_sphere, equator_ring, axis_line)
        return earth_group

    def position_earth_at_angle(
        self, earth_group, orbit_center, angle, with_tilt=False
    ):
        """Position Earth at a specific orbital angle"""
        earth_group.move_to(
            orbit_center
            + self.ORBIT_RADIUS * np.array([np.cos(angle), np.sin(angle), 0])
        )

    def construct(self):
        self.setup_constants()
        orbit_center = ORIGIN

        # Camera angle (Fixed throughout)
        self.set_camera_orientation(phi=self.PHI, theta=self.THETA)

        # Create orbital system
        orbital_system = self.create_orbital_system(orbit_center)

        # ===== PHASE 1: NO TILT =====
        title1 = Text("Phase 1: No Axial Tilt", font_size=32).to_edge(UP)
        self.add_fixed_in_frame_mobjects(title1)
        self.play(Write(title1), Create(orbital_system), run_time=2)

        # Create Earth No Tilt
        earth_no_tilt = self.create_earth(with_tilt=False)
        # Position at start (angle 0)
        self.position_earth_at_angle(earth_no_tilt, orbit_center, 0)
        self.add(earth_no_tilt)
        self.play(FadeIn(earth_no_tilt))

        # Full Orbit
        orbit_tracker = ValueTracker(0)

        def orbit_updater_no_tilt(mob):
            angle = orbit_tracker.get_value()
            self.position_earth_at_angle(mob, orbit_center, angle, with_tilt=False)

        earth_no_tilt.add_updater(orbit_updater_no_tilt)

        # Label: "No axial tilt - axis perpendicular to orbit"
        label1 = Text("Axis perpendicular to orbit", font_size=24).to_edge(DOWN)
        self.add_fixed_in_frame_mobjects(label1)
        self.play(Write(label1))

        self.play(orbit_tracker.animate.set_value(2 * PI), run_time=6, rate_func=linear)
        earth_no_tilt.remove_updater(orbit_updater_no_tilt)
        self.wait(0.5)

        # Transition to Phase 2
        self.play(FadeOut(title1), FadeOut(label1), FadeOut(earth_no_tilt), run_time=1)

        # ===== PHASE 2: WITH TILT =====
        title2 = Text("Phase 2: With 23.5° Axial Tilt", font_size=32).to_edge(UP)
        self.add_fixed_in_frame_mobjects(title2)
        self.play(Write(title2))

        # Create Earth With Tilt
        earth_with_tilt = self.create_earth(with_tilt=True)
        # Start from 0
        orbit_tracker.set_value(0)
        self.position_earth_at_angle(earth_with_tilt, orbit_center, 0, with_tilt=True)
        self.add(earth_with_tilt)
        self.play(FadeIn(earth_with_tilt))

        def orbit_updater_with_tilt(mob):
            angle = orbit_tracker.get_value()
            self.position_earth_at_angle(mob, orbit_center, angle, with_tilt=True)

        earth_with_tilt.add_updater(orbit_updater_with_tilt)

        # Label: "With 23.5° axial tilt" (Already in title)

        self.play(orbit_tracker.animate.set_value(2 * PI), run_time=6, rate_func=linear)
        earth_with_tilt.remove_updater(orbit_updater_with_tilt)
        self.wait(0.5)

        # ===== PHASE 3: TRANSFORM TO TILTED ORBITAL PLANE =====
        self.play(FadeOut(title2))
        title3 = Text("Phase 3: Tilted Orbital Plane", font_size=32).to_edge(UP)
        self.add_fixed_in_frame_mobjects(title3)
        self.play(Write(title3))

        # Create sugar group for system (Sun + Orbit)
        # We separate Earth because we want to animate its rotation independently/differently
        system_group = VGroup(orbital_system, earth_with_tilt)

        # Smooth continuous transition
        # "Rotate the entire orbital system (sun, orbit circle, Earth) around X-axis by 23.44°"
        # "Also rotate Earth itself so its axis appears flat/horizontal relative to the camera"

        # 1. Separate Earth from system temporarily so we can master its rotation
        system_group.remove(earth_with_tilt)
        orbit_system_only = system_group

        # 2. Create the TARGET Earth (Horizontal Axis)
        # Horizontal Axis = RIGHT (X-axis).
        # We create a fresh Earth, rotate its axis to be Horizontal.
        # Use create_earth helper.
        # Base Earth has Axis UP. with_tilt=False.
        earth_target = self.create_earth(with_tilt=False)
        # Rotate -90 around Z (OUT) to make Axis (UP) become Horizontal (RIGHT)
        # Point (0, 1, 0) -> (1, 0, 0)
        earth_target.rotate(-90 * DEGREES, axis=OUT)
        # Move to the correct position (same as current earth: R, 0, 0)
        earth_target.move_to(orbit_center + np.array([self.ORBIT_RADIUS, 0, 0]))

        # 3. Animate
        # System rotates TILT around RIGHT.
        # Earth transforms to Target (Horizontal).
        # Note: EARTH POSITION stays at (R, 0, 0) during X-axis rotation of system!
        # So we don't need to move earth_target, just rotate it.

        self.play(
            Rotate(
                orbit_system_only, angle=self.TILT, axis=RIGHT, about_point=orbit_center
            ),
            Transform(earth_with_tilt, earth_target),
            run_time=3,
        )

        # 4. Independent Update loop
        self.add(earth_with_tilt)  # Ensure clean state

        tracker_phase3 = ValueTracker(0)

        def orbit_updater_3d(mob):
            angle = tracker_phase3.get_value()

            # Position on Circle in XY plane (Original Orbit)
            x = self.ORBIT_RADIUS * np.cos(angle)
            y = self.ORBIT_RADIUS * np.sin(angle)
            z = 0

            # Apply the System Tilt (Rotation around X by TILT)
            t = self.TILT
            x_new = x
            y_new = y * np.cos(t) - z * np.sin(t)
            z_new = y * np.sin(t) + z * np.cos(t)

            mob.move_to(orbit_center + np.array([x_new, y_new, z_new]))
            # Orientation is fixed (Horizontal) by the Transform so no rotation needed per frame.

        earth_with_tilt.add_updater(orbit_updater_3d)

        note3 = Text(
            "Earth axis is horizontal; Orbit is tilted 23°", font_size=24, color=YELLOW
        ).to_edge(DOWN)
        self.add_fixed_in_frame_mobjects(note3)
        self.play(Write(note3))

        self.play(
            tracker_phase3.animate.set_value(2 * PI), run_time=8, rate_func=linear
        )

        earth_with_tilt.remove_updater(orbit_updater_3d)
        self.wait(2)

        self.play(
            FadeOut(orbit_system_only),
            FadeOut(earth_with_tilt),
            FadeOut(title3),
            FadeOut(note3),
            run_time=2,
        )
