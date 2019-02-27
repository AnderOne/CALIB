%Выполняет анимацию результатов моделирования, представленных в виде сеточных полей:
%fname -- имя файла с контурами; min_t -- номер начального кадра; max_t -- номер последнего кадра; jmp_t -- число пропускаемых кадров;
function draw_grid(fname, min_t, max_t, jmp_t, gif)

	if (nargin < 1) error('Input file is not specified!'); end
	if (nargin < 2) min_t = 0; end
	if (nargin < 3) max_t = min_t; end
	if (nargin < 4) jmp_t = 1; end
	if (nargin < 5) gif = 0; end

	fid = fopen(fname, 'rb');
	%Считываем заголовок:
	min_x = fread(fid, 1, 'double'); min_y = fread(fid, 1, 'double'); max_x = fread(fid, 1, 'double'); max_y = fread(fid, 1, 'double');
	num_t = fread(fid, 1, 'uint32'); num_t = min(num_t, max_t);
	num_x = fread(fid, 1, 'uint32'); len_x = (max_x - min_x) / (num_x - 1);
	num_y = fread(fid, 1, 'uint32'); len_y = (max_y - min_y) / (num_y - 1);
	for t = 0 : num_t
		%Читаем очередной кадр из файла:
		if ((t >= min_t) && ~mod(t, jmp_t))
			val_z = fread(fid, [num_x, num_y], 'double');
		else
			%Пропускаем лишние кадры:
			fseek(fid, 8 * num_x * num_y, 'cof');
			continue
		end
		%Выполняем отрисовку поля:
if (true)
		pcolor(min_x + (0 : num_x - 1) * len_x, min_y + (0 : num_y - 1) * len_y, val_z.');
		colormap('default')
else
		surfl(min_x + (0 : num_x - 1) * len_x, min_y + (0 : num_y - 1) * len_y, val_z.');
		colormap('gray')
end
		shading flat
		title(['timestep: ', num2str(t)]);
		axis([min_x, max_x, min_y, max_y]);
		mov = getframe();
		%drawnow
		pause(0.1);
		hold off
		%...
		if (~gif) continue; end
		[im, map] = rgb2ind(mov.cdata, 256);
		if (t > min_t)
			imwrite(im, map, 'grid.gif', 'DelayTime', 0.1, 'WriteMode', 'Append');
		else
			%imwrite(im, map, 'grid.gif', 'DelayTime', 0.1);%, 'LoopCount', inf);
			imwrite(im, map, 'grid.gif', 'DelayTime', 0.1, 'LoopCount', inf);
		end
	end
	%if (gif) imwrite(im, map, 'grid.png'); end

	fclose(fid);
end
