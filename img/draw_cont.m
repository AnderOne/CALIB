%Выполняет анимацию результатов моделирования, представленных в виде карты изолиний:
%fname -- имя файла с контурами; min_t -- номер начального кадра; max_t -- номер последнего кадра; jmp_t -- число пропускаемых кадров;
function draw_cont(fname, min_t, max_t, jmp_t, rec, gif)

	if (nargin < 1) error('Input file is not specified!'); end
	if (nargin < 2) min_t = 0; end
	if (nargin < 3) max_t = min_t; end
	if (nargin < 4) jmp_t = 1; end
	if (nargin < 5) rec = []; end
	if (nargin < 6) gif = 0; end
	if ((length(rec) ~= 0) && (length(rec) ~= 4))
		error('Incorrect frame parameter!');
	end

	fid = fopen(fname, 'rb');
	%Считываем заголовок:
	min_x = fread(fid, 1, 'double'); min_y = fread(fid, 1, 'double'); max_x = fread(fid, 1, 'double'); max_y = fread(fid, 1, 'double');
	num_t = fread(fid, 1, 'uint32'); num_t = min(num_t, max_t);
	hal_x = (max_x - min_x) / 2; hal_y = (max_y - min_y) / 2;
	if (length(rec) == 0)
		rec = [min_x, max_x, min_y, max_y];
	end
	for t = 0 : num_t
		if ((t >= min_t) && ~mod(t, jmp_t))
			%Читаем очередной кадр из файла:
			num_p = fread(fid, 1, 'uint32');
			NODEX = fread(fid, num_p, 'double'); NODEY = fread(fid, num_p, 'double');
			num_c = fread(fid, 1, 'uint32');
			CONTI = fread(fid, num_c, 'uint32'); CONTN = fread(fid, num_c, 'uint32'); CONTS = fread(fid, num_c, 'uchar');
			num_s = fread(fid, 1, 'uint32');
			STEPI = fread(fid, num_s, 'uint32'); STEPN = fread(fid, num_s, 'uint32');
		else
			%Пропускаем лишние кадры:
			num_p = fread(fid, 1, 'uint32'); fseek(fid, 2 * 8 * num_p, 'cof');
			num_c = fread(fid, 1, 'uint32'); fseek(fid, 2 * 4 * num_c, 'cof');  fseek(fid, num_c, 'cof');
			num_s = fread(fid, 1, 'uint32'); fseek(fid, 2 * 4 * num_s, 'cof');
			continue
		end
		%Выполняем отрисовку контуров:
		RGB = [0, 0, 0]; dc = 0.5 * [0, 1, 2]; hc = 1.0 / num_s;
		for l = 1 : num_s
			c1 = STEPI(l) + 1; c2 = (c1 + STEPN(l)) - 1;
			for i = c1 : c2
				k1 = CONTI(i) + 1; k2 = (k1 + CONTN(i)) - 1;
				i1 = k1; f = false;
				for k = k1 : (k2 - 1)
					x1 = NODEX(k); x2 = NODEX(k + 1);
					dx = x2 - x1; if (abs(dx) >= hal_x); f = true; dx = dx - 2.0 * hal_x * fix(dx / hal_x); end
					y1 = NODEY(k); y2 = NODEY(k + 1);
					dy = y2 - y1; if (abs(dy) >= hal_y); f = true; dy = dy - 2.0 * hal_y * fix(dy / hal_y); end
					if (~f); continue; end
					plot(NODEX(i1 : k), NODEY(i1 : k), 'Color', RGB);%, 'Marker', '.');
					hold on
					plot([x1, x1 + dx], [y1, y1 + dy], 'Color', RGB);%, 'Marker', '.');
					plot([x2, x2 - dx], [y2, y2 - dy], 'Color', RGB);%, 'Marker', '.');
					plot([x1, x1 + dx], [y2 - dy, y2], 'Color', RGB);%, 'Marker', '.');
					plot([x2, x2 - dx], [y1 + dy, y1], 'Color', RGB);%, 'Marker', '.');
					i1 = k + 1;
					f = false;
				end
				plot(NODEX(i1 : k2), NODEY(i1 : k2), 'Color', RGB);%, 'Marker', '.');
				hold on
				if (CONTS(i))
					x1 = NODEX(k2); x2 = NODEX(k1);
					dx = x2 - x1; if (abs(dx) >= hal_x); dx = dx - 2.0 * hal_x * fix(dx / hal_x); end
					y1 = NODEY(k2); y2 = NODEY(k1);
					dy = y2 - y1; if (abs(dy) >= hal_y); dy = dy - 2.0 * hal_y * fix(dy / hal_y); end
					plot([x1, x1 + dx], [y1, y1 + dy], 'Color', RGB);%, 'Marker', '.');
					plot([x2, x2 - dx], [y2, y2 - dy], 'Color', RGB);%, 'Marker', '.');
					plot([x1, x1 + dx], [y2 - dy, y2], 'Color', RGB);%, 'Marker', '.');
					plot([x2, x2 - dx], [y1 + dy, y1], 'Color', RGB);%, 'Marker', '.');
					i1 = k + 1;
					f = false;
				end
				hold on
			end
			RGB = 2 * (l * hc >= dc).* (l * hc - dc);
			RGB = (RGB <= 1).* RGB + (RGB > 1);
		end
		title(['timestep: ', num2str(t), '; nk(', num2str(num_p), '); nc (', num2str(num_c), ')']);
		axis(rec);
		mov = getframe();
		%drawnow
		pause(0.1);
		hold off
		%...
		if (~gif) continue; end
		[im, map] = rgb2ind(mov.cdata, 256);
		if (t > min_t)
			imwrite(im, map, 'cont.gif', 'DelayTime', 0.1, 'WriteMode', 'Append');
		else
			imwrite(im, map, 'cont.gif', 'DelayTime', 0.1);%, 'LoopCount', inf);
			%imwrite(im, map, 'cont.gif', 'DelayTime', 0.1, 'LoopCount', inf);
		end
	end
	%if (gif) imwrite(im, map, 'cont.png'); end

	fclose(fid);
end
