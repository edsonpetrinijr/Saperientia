from manimlib import *
import os
from pathlib import Path

class VideoMobject(Group):
    def __init__(
        self,
        directory,
        frame_rate=30,
        loop=False,
        start_frame=0,
        end_frame=None,
        **kwargs
    ):
        """
        Create a video-like object from a directory of image frames.
        
        Parameters:
        -----------
        directory : str
            Path to directory containing image frames (relative or absolute)
        frame_rate : int
            Frames per second for playback (default: 30)
        loop : bool
            Whether to loop the video (default: False)
        start_frame : int
            Frame index to start from (default: 0)
        end_frame : int or None
            Frame index to end at, None for all frames (default: None)
        """
        super().__init__(**kwargs)
        
        self.directory = directory
        self.frame_rate = frame_rate
        self.loop = loop
        self.start_frame = start_frame
        self.end_frame = end_frame
        
        # Load all image paths
        self.frame_paths = self._load_frame_paths()
        
        if not self.frame_paths:
            raise ValueError(f"No image files found in directory: {directory}")
        
        # Apply frame range
        if end_frame is not None:
            self.frame_paths = self.frame_paths[start_frame:end_frame]
        else:
            self.frame_paths = self.frame_paths[start_frame:]
        
        self.total_frames = len(self.frame_paths)
        self.current_frame_index = 0
        self.time_counter = 0
        
        # Create initial image
        self.current_image = ImageMobject(self.frame_paths[0])
        self.add(self.current_image)
        
        # Add updater to animate frames
        self.add_updater(self._update_frame)
        
    def _load_frame_paths(self):
        """Load and sort all image file paths from the directory."""
        path = Path(self.directory)
        
        if not path.exists():
            raise FileNotFoundError(f"Directory not found: {self.directory}")
        
        # Common image extensions
        extensions = ('.png', '.jpg', '.jpeg', '.bmp', '.gif', '.tiff')
        
        # Get all image files
        frame_files = [
            f for f in path.iterdir() 
            if f.is_file() and f.suffix.lower() in extensions
        ]
        
        # Sort files naturally (handles frame_001.png, frame_002.png, etc.)
        frame_files.sort(key=lambda x: x.name)
        
        return [str(f) for f in frame_files]
    
    def _update_frame(self, mob, dt):
        """Update the displayed frame based on elapsed time."""
        self.time_counter += dt
        
        # Calculate which frame should be displayed
        target_frame = int(self.time_counter * self.frame_rate)
        
        # Handle looping
        if self.loop:
            target_frame = target_frame % self.total_frames
        else:
            target_frame = min(target_frame, self.total_frames - 1)
        
        # Update the image if frame changed
        if target_frame != self.current_frame_index:
            self._set_frame(target_frame)
            self.current_frame_index = target_frame
    
    def _set_frame(self, frame_idx):
        """Set the displayed image to a specific frame."""
        if 0 <= frame_idx < self.total_frames:
            # Store current properties
            height = self.current_image.get_height()
            width = self.current_image.get_width()
            center = self.current_image.get_center()
            
            # Remove old image
            self.remove(self.current_image)
            
            # Create and add new image
            self.current_image = ImageMobject(self.frame_paths[frame_idx])
            self.current_image.set_height(height)
            self.current_image.move_to(center)
            self.add(self.current_image)
    
    def reset(self):
        """Reset the video to the first frame."""
        self.current_frame_index = 0
        self.time_counter = 0
        self._set_frame(0)
        
    def seek(self, frame_index):
        """Jump to a specific frame."""
        self.current_frame_index = frame_index
        self.time_counter = frame_index / self.frame_rate
        self._set_frame(int(frame_index))