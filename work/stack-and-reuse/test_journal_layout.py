"""Regression checks for panel lettering without changing plotted coordinates."""
import unittest
import matplotlib.pyplot as plt
from galstyle import journal_layout, panel


class JournalLayoutTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_unlettered_plot_gets_clear_label_without_changing_limits(self):
        fig, ax = plt.subplots(figsize=(7.2, 3))
        ax.plot([2013, 2025], [0, 25])
        ax.set_yticks([0, 25], ['0%', '25%'])
        limits = ax.get_xlim(), ax.get_ylim(), ax.get_position().bounds
        bounds = journal_layout(fig)
        label = ax.texts[-1]
        self.assertEqual(label.get_text(), 'a')
        self.assertEqual(limits, (ax.get_xlim(), ax.get_ylim(), ax.get_position().bounds))
        renderer = fig.canvas.get_renderer()
        self.assertGreater(label.get_window_extent(renderer).y0, ax.bbox.y1)
        self.assertLessEqual(bounds.x0, 0)
        self.assertGreaterEqual(bounds.x1, fig.get_figwidth())

    def test_row_letters_align_with_different_title_heights(self):
        fig = plt.figure(figsize=(7.2, 3))
        a = fig.add_axes([.1, .15, .35, .6])
        b = fig.add_axes([.65, .15, .3, .6])
        a.set_title('One line')
        b.set_title('Two\nlines')
        panel(a, 'a', y=1.02)
        panel(b, 'b', y=1.12)
        journal_layout(fig)
        self.assertAlmostEqual(a.texts[-1].get_position()[1], b.texts[-1].get_position()[1])
        renderer = fig.canvas.get_renderer()
        self.assertFalse(b.texts[-1].get_window_extent(renderer).overlaps(b.title.get_window_extent(renderer)))

    def test_lower_right_panel_does_not_jump_to_left_gutter(self):
        fig = plt.figure(figsize=(7.2, 5))
        a = fig.add_axes([.1, .1, .35, .75])
        b = fig.add_axes([.7, .55, .25, .3])
        c = fig.add_axes([.7, .1, .25, .3])
        for ax, letter in zip([a, b, c], 'abc'):
            panel(ax, letter)
        journal_layout(fig)
        self.assertGreater(c.texts[-1].get_position()[0], .5)
        self.assertAlmostEqual(b.texts[-1].get_position()[0], c.texts[-1].get_position()[0])


if __name__ == '__main__':
    unittest.main()
