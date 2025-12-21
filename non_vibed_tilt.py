from manim import *


class SiderealVsSolarNoTilt(Scene):
    def setup(self):
        # Constants
        self.EARTH_RADIUS = 0.5
        self.ROTATION_DURATION = 5
        self.NUM_SPINS = 1.0
        self.SOLAR_SYSTEM_RADIUS = 3.0
        self.ORBIT_DURATION = 120.0  # Arbitrary slow orbit for now
        self.POINTER_LENGTH = 0.7

        self.sun = Circle(radius=1).set_color(YELLOW).set_opacity(0.3)

        self.planetary_rim = (
            Circle(radius=self.SOLAR_SYSTEM_RADIUS).set_color(GRAY).set_stroke(width=1)
        )
        self.planetary_rim.set_fill(GRAY, opacity=0.2)

        self.earth = (
            Circle(radius=self.EARTH_RADIUS)
            .set_color(BLUE_E)
            .set_stroke(width=1)
            .set_opacity(0.5)
        )
        earth_center = np.array([0, self.SOLAR_SYSTEM_RADIUS, 0])
        self.earth.move_to(earth_center)

        # Line starts from the edge of the earth (center - radius in y)
        line_start = earth_center - np.array([0, self.EARTH_RADIUS, 0])
        # Original line was ~0.7 units long (2.5 to 1.8)
        line_end = line_start - np.array([0, self.POINTER_LENGTH, 0])
        self.earth_line = Line(line_start, line_end)
        self.earth_line.set_color(YELLOW_A)

        self.earth_and_stick = VGroup(self.earth, self.earth_line)

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
        self.earth_and_stick.save_state()
        # Initialize attribute after save_state so it is not part of the 'saved' state?
        # Actually simplest to just set it.
        self.earth_and_stick.old_rotation_angle = 0

        self.play(Create(self.sun))
        self.play(
            Create(self.planetary_rim),
            Create(self.earth_and_stick),
        )
        self.add(self.spin_label, self.spin_value)

        def run_day_cycle(
            target_angle: float,
            run_time: float = self.ROTATION_DURATION,
            show_solar_arrow: bool = False,
        ) -> None:
            # Reset
            self.earth_and_stick.restore()
            self.orbit_tracker.set_value(PI / 2)
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
