import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.image as mpimg
import glob
import re

# ============================================================
# 1. Collect JPG files named like:
#    rho_xy_it_0.0.jpg
#    rho_xy_it_1024.0.jpg
#    rho_xy_it_2048.0.jpg
# ============================================================
file_pattern = "rho_xy_it_*.jpg"

files = sorted(
    glob.glob(file_pattern),
    key=lambda name: float(re.search(r"it_([0-9.]+)\.jpg$", name).group(1))
)

if len(files) == 0:
    raise RuntimeError(f"No JPG files found matching pattern {file_pattern}")


# ============================================================
# 2. Extract iteration numbers for the title
# ============================================================
def extract_iteration(fname):
    m = re.search(r"it_([0-9.]+)\.jpg$", fname)
    return float(m.group(1)) if m else None

iterations = [extract_iteration(f) for f in files]

MSUN_TO_MS = 4.925490947e-3   # 1 M_sun in ms
DT_MSU = 0.0625
# 1024 iter =  64.0 Msun                 # timestep in M_sun
tmerg = 12.0

times_ms = [it * DT_MSU * MSUN_TO_MS - tmerg for it in iterations]

#iterations = [extract_iteration(f) for f in files]


# ============================================================
# 3. Load first image and setup figure
# ============================================================
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

fig, ax = plt.subplots(figsize=(6,6))

img0 = mpimg.imread(files[0])
im = ax.imshow(img0)
ax.axis("off")

tmerg = 12.0

title = ax.set_title(r"$\rho,\ t = %.3f\ \mathrm{ms}$" % times_ms[0])



# ============================================================
# 4. Update function for animation
# ============================================================
def update(frame):
    img = mpimg.imread(files[frame])
    im.set_data(img)

    #it = iterations[frame]
    #title.set_text(r"$\rho,\ \mathrm{iteration}\ %d$" % it)
    t_ms = times_ms[frame]
    title.set_text(r"$\rho,\ t-t_{\rm mer} = %.3f\ \mathrm{ms}$" % t_ms)

    return im, title


# ============================================================
# 5. Animate
# ============================================================
anim = FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=160, #120,
    blit=True
)

# ============================================================
# 6. Save animation
# ============================================================
anim.save("rho_animation.mp4", writer="ffmpeg", dpi=150)
# anim.save("rho_animation.gif", writer="pillow", dpi=150)

# plt.show()

