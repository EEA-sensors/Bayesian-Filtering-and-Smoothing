import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_pendulum", "plot_car_trajectory"]


def plot_pendulum(timeline, y, x1, label1, x2=None, label2=None):
    fig, axes = plt.subplots(ncols=2, figsize=(22, 10))
    axes[1].scatter(timeline, y, marker="o", label="Measurements", color="red", alpha=0.66)
    axes[1].plot(timeline, np.sin(x1[:, 0]), linestyle="dashdot", label=label1, color="blue")

    axes[0].plot(x1[:, 0], x1[:, 1], label=label1, color="blue")

    if x2 is not None:
        axes[1].plot(timeline, np.sin(x2[:, 0]), linestyle="dashdot", label=label2, color="orange")
        axes[0].plot(x2[:, 0], x2[:, 1], label=label2, color="orange")
    else:
        axes[0].scatter(x1[0, 0], x1[0, 1], marker="x", color="orange", s=500)

    axes[0].set_xlabel("$x_0(t)$")
    axes[0].set_ylabel("$x_1(t)$")

    axes[1].set_xlabel("$t$")
    axes[1].set_ylabel("$\sin(x_0(t))$")

    axes[0].legend()
    axes[1].legend()


def plot_car_trajectory(y, x1, label1, x2=None, label2=None):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(y[:, 0], y[:, 1], marker="o", label="Measurements", color="red")
    ax.plot(x1[:, 0], x1[:, 1], label=label1, color="blue")
    if x2 is None:
        ax.scatter(x1[0, 0], x1[0, 1], marker="x", color="orange", s=500)
    else:
        ax.plot(x2[:, 0], x2[:, 1], label=label2, color="orange")
    _ = ax.legend()
    _ = ax.set_xlabel("${\it x}_1$")
    _ = ax.set_ylabel("${\it x}_2$")
