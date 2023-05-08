from manim import *
import numpy as np
from scipy import integrate
import sympy


def f(x):
    return np.cos(2*x)


class Center(Scene):
    def construct(self):
        ax = Axes(
            x_range=[-0.1, 1.7], y_range=[-1.5, 1.5], tips=True,
            axis_config={
                "include_numbers": True,
                "color": WHITE
            },
        )
        labels = ax.get_axis_labels(x_label="x", y_label="y=cos(2x)")
        self.play(Create(ax), Create(labels))
        self.add(ax, labels)

        graph = ax.plot(lambda x: np.cos(2*x), x_range=[0, np.pi/2])
        graph.show()
        self.play(Create(graph))

        riemann_areas = [
            ax.get_riemann_rectangles(graph, x_range=[0, np.pi/4], input_sample_type="center", color=('#ff69b4', '#ffff00'),
                                      dx= np.pi / (4*i), fill_opacity=0.8) for i in range(2, 110, 2)]
        tex = [Tex(f"$n={i}$").to_corner(RIGHT + UP) for i in range(2, 110, 2)]
        tex_last = Tex(r"$n \to \infty$").to_corner(RIGHT + UP)
        type = [Tex(r'Center').next_to(
            tex[(n - 2) // 2], DOWN) for n in range(2, 110, 2)]

        real_area, e = integrate.quad(f, 0.0, np.pi/4)

        current_area = [Tex(r'$\sigma = $' + "{0:.6f}".format(
            (sum([f((np.pi * i / (4*n) + np.pi * (i - 1) / (4*n)) / 2) * np.pi / (4*n) for i in range(1, n + 1)])))).to_edge(UP) for n in
                        range(2, 110, 2)]

        current_R_n = [Tex(r'$Real \left|R_n \right| = $' + "{0:.6f}".format(
            abs(real_area - sum([f((np.pi * i / (4 * n) + np.pi * (i - 1) / (4 * n)) / 2) * np.pi / (4 * n) for i in
                                  range(1, n + 1)])))).next_to(
            current_area[(n - 2)//2], DOWN) for n in range(2, 110, 2)]
        theory_R_n = [Tex(r'$Theory \left|R_n \right| = $' + "{0:.6f}".format((np.pi ** 3) / (384 * (n ** 2)))).next_to(current_R_n[(n - 2)//2], DOWN) for n in
                      range(2, 110, 2)]

        self.play(Create(riemann_areas[0]))
        self.add(tex[0])
        self.add(type[0])
        self.add(current_area[0])
        self.add(current_R_n[0])
        self.add(theory_R_n[0])

        self.wait(0.8)

        for i in range(0, len(riemann_areas), 4):
            self.play(Transform(riemann_areas[0], riemann_areas[i]), Transform(tex[0], tex[i]),
                      Transform(current_area[0], current_area[i]), Transform(current_R_n[0], current_R_n[i]),
                      Transform(theory_R_n[0], theory_R_n[i]), run_time=0.2)
            self.wait(0.8)

        area = ax.get_area(graph, x_range=[0, np.pi /4], opacity=1, color=('#ff69b4', '#ffff00'))
        self.remove(current_R_n[0], theory_R_n[0])
        self.play(Transform(riemann_areas[0], area), Transform(tex[0], tex_last),
                  Transform(current_area[0], Tex(r'$\int\limits_0^{\frac{\pi}{4}} cos(2x)\, dx = $' + str(round(real_area, 1))).to_edge(UP)))
        self.wait(2)


class Right(Scene):
    def construct(self):
        ax = Axes(
            x_range=[-0.1, 1.7], y_range=[-1.5, 1.5], tips=True,
            axis_config={
                "include_numbers": True,
                "color": WHITE
            },
        )
        labels = ax.get_axis_labels(x_label="x", y_label="y=cos(2x)")
        self.play(Create(ax), Create(labels))
        self.add(ax, labels)

        graph = ax.plot(lambda x: np.cos(2*x), x_range=[0, np.pi/2])
        graph.show()
        self.play(Create(graph))

        riemann_areas = [
            ax.get_riemann_rectangles(graph, x_range=[0, np.pi/4], input_sample_type="right", color=('#fefe22', '#ff4500'),
                                      dx= np.pi / (4*i), fill_opacity=0.8) for i in range(2, 110, 2)]
        tex = [Tex(f"$n={i}$").to_corner(RIGHT + UP) for i in range(2, 110, 2)]
        tex_last = Tex(r"$n \to \infty$").to_corner(RIGHT + UP)
        type = [Tex(r'Right').next_to(
            tex[(n - 2) // 2], DOWN) for n in range(2, 110, 2)]

        real_area, e = integrate.quad(f, 0.0, np.pi/4)

        current_area = [Tex(r'$\sigma = $' + "{0:.6f}".format(
            (sum([f(np.pi * i / (4*n)) * np.pi / (4*n) for i in range(1, n + 1)])))).to_edge(UP) for n in
                        range(2, 110, 2)]

        current_R_n = [Tex(r'$Real \left|R_n \right| = $' + "{0:.6f}".format(
            abs(real_area - sum([f(np.pi * i / (4 * n)) * np.pi / (4 * n) for i in
                                  range(1, n + 1)])))).next_to(
            current_area[(n - 2)//2], DOWN) for n in range(2, 110, 2)]
        theory_R_n = [Tex(r'$Theory \left|R_n \right| = $' + "{0:.6f}".format((np.pi ** 2) / (16 * n))).next_to(current_R_n[(n - 2)//2], DOWN) for n in
                      range(2, 110, 2)]

        self.play(Create(riemann_areas[0]))
        self.add(tex[0])
        self.add(type[0])
        self.add(current_area[0])
        self.add(current_R_n[0])
        self.add(theory_R_n[0])

        self.wait(0.8)

        for i in range(0, len(riemann_areas), 4):
            self.play(Transform(riemann_areas[0], riemann_areas[i]), Transform(tex[0], tex[i]),
                      Transform(current_area[0], current_area[i]), Transform(current_R_n[0], current_R_n[i]),
                      Transform(theory_R_n[0], theory_R_n[i]), run_time=0.2)
            self.wait(0.8)

        area = ax.get_area(graph, x_range=[0, np.pi /4], opacity=1, color=('#fefe22', '#ff4500'))
        self.remove(current_R_n[0], theory_R_n[0])
        self.play(Transform(riemann_areas[0], area), Transform(tex[0], tex_last),
                  Transform(current_area[0], Tex(r'$\int\limits_0^{\frac{\pi}{4}} cos(2x)\, dx = $' + str(round(real_area, 1))).to_edge(UP)))
        self.wait(2)


class Left(Scene):
    def construct(self):
        ax = Axes(
            x_range=[-0.1, 1.7], y_range=[-1.5, 1.5], tips=True,
            axis_config={
                "include_numbers": True,
                "color": WHITE
            },
        )
        labels = ax.get_axis_labels(x_label="x", y_label="y=cos(2x)")
        self.play(Create(ax), Create(labels))
        self.add(ax, labels)

        graph = ax.plot(lambda x: np.cos(2*x), x_range=[0, np.pi/2])
        graph.show()
        self.play(Create(graph))

        riemann_areas = [
            ax.get_riemann_rectangles(graph, x_range=[0, np.pi/4], input_sample_type="left", color=('#add8e6', '#8b008b'),
                                      dx= np.pi / (4*i), fill_opacity=0.8) for i in range(2, 110, 2)]
        tex = [Tex(f"$n={i}$").to_corner(RIGHT + UP) for i in range(2, 110, 2)]
        tex_last = Tex(r"$n \to \infty$").to_corner(RIGHT + UP)
        type = [Tex(r'Left').next_to(
            tex[(n - 2)//2], DOWN) for n in range(2, 110, 2)]
        real_area, e = integrate.quad(f, 0.0, np.pi/4)

        current_area = [Tex(r'\qquad $\sigma = $' + "{0:.6f}".format(
            (sum([f( np.pi * (i - 1) / (4*n)) * np.pi / (4*n) for i in range(1, n + 1)])))).to_edge(UP) for n in
                        range(2, 110, 2)]

        current_R_n = [Tex(r'\qquad $Real \left|R_n \right| = $' + "{0:.6f}".format(
            abs(real_area - sum([f(np.pi * (i - 1) / (4 * n)) * np.pi / (4 * n) for i in
                                  range(1, n + 1)])))).next_to(
            current_area[(n - 2)//2], DOWN) for n in range(2, 110, 2)]
        theory_R_n = [Tex(r'\qquad $Theory \left|R_n \right| = $' + "{0:.6f}".format((np.pi ** 2) / (16 * n))).next_to(current_R_n[(n - 2)//2], DOWN) for n in
                      range(2, 110, 2)]

        self.play(Create(riemann_areas[0]))
        self.add(tex[0])
        self.add(type[0])
        self.add(current_area[0])
        self.add(current_R_n[0])
        self.add(theory_R_n[0])

        self.wait(0.8)

        for i in range(0, len(riemann_areas), 4):
            self.play(Transform(riemann_areas[0], riemann_areas[i]), Transform(tex[0], tex[i]),
                      Transform(current_area[0], current_area[i]), Transform(current_R_n[0], current_R_n[i]),
                      Transform(theory_R_n[0], theory_R_n[i]), run_time=0.2)
            self.wait(0.8)

        area = ax.get_area(graph, x_range=[0, np.pi /4], opacity=1, color=('#add8e6', '#8b008b'))
        self.remove(current_R_n[0], theory_R_n[0])
        self.play(Transform(riemann_areas[0], area), Transform(tex[0], tex_last),
                  Transform(current_area[0], Tex(r'$\int\limits_0^{\frac{\pi}{4}} cos(2x)\, dx = $' + str(round(real_area, 1))).to_edge(UP)))
        self.wait(2)
